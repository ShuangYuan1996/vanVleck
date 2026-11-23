"""Main window: file loading, dataset switching, display options, coordination-sphere panel."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor, QDoubleValidator, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDockWidget, QFileDialog, QMenu,
                               QTreeWidget, QTreeWidgetItem,
                               QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
                               QTabWidget, QTextBrowser, QToolBar, QVBoxLayout, QWidget)

from ..coordination import CoordinationSphere, sphere_for_clicked_atom
from .. import __app_name__, __url__, __version__
from .. import symmetry as sym
from .. import vanvleck as vv
from .. import distortion as dist
from ..structure import Structure, load_file
from .compare import CompareDialog
from .masses import MassEditor
from .qchart import QChartDialog, export_q_dialog
from .viewer import STYLES, MoleculeView

FILE_FILTER = "Structure files (*.cif *.xyz);;CIF files (*.cif);;XYZ files (*.xyz);;All files (*)"


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("vanVleck — Spin-Crossover Structure Viewer")
        self.resize(1300, 820)
        self.setAcceptDrops(True)

        self.structures: list[Structure] = []      # datasets of all files in the workspace (flat list)
        self.files: list[dict] = []                # [{'path': path, 'structures': [...]}]
        self.current: Structure | None = None
        self.current_sphere: CoordinationSphere | None = None   # used by later analysis steps
        self.current_symmetry: sym.SymmetryResult | None = None
        self.current_vanvleck: vv.VanVleckResult | None = None
        self.qdialog: QChartDialog | None = None
        self.compare_dialog: CompareDialog | None = None

        self.view = MoleculeView()
        self.view.atomClicked.connect(self._on_atom_clicked)
        self.view.atomHovered.connect(self._on_atom_hovered)

        self._build_toolbar()
        self._build_menu()
        self._build_side_panel()
        self._build_navigator()

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.view)
        splitter.addWidget(self.side)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([880, 420])
        self.setCentralWidget(splitter)
        self.statusBar().showMessage("Ready — open a .cif or .xyz file")

    # ------------------------------------------------------------------ UI
    def _build_toolbar(self):
        tb = QToolBar("Main toolbar")
        tb.setMovable(False)
        self.addToolBar(tb)

        self.act_open = QAction("Add files…", self)
        self.act_open.setShortcut(QKeySequence.Open)
        self.act_open.triggered.connect(self.open_dialog)
        tb.addAction(self.act_open)
        tb.addSeparator()

        tb.addWidget(QLabel(" Dataset: "))
        self.combo_data = QComboBox()
        self.combo_data.setMinimumWidth(380)
        self.combo_data.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.combo_data.currentIndexChanged.connect(self._on_dataset_changed)
        tb.addWidget(self.combo_data)
        tb.addSeparator()

        tb.addWidget(QLabel(" Style: "))
        self.combo_style = QComboBox()
        for key, name in STYLES.items():
            self.combo_style.addItem(name, key)
        self.combo_style.currentIndexChanged.connect(
            lambda: self.view.set_style(self.combo_style.currentData()))
        tb.addWidget(self.combo_style)
        tb.addSeparator()

        def check(text, attr, default):
            cb = QCheckBox(text)
            cb.setChecked(default)
            cb.toggled.connect(lambda v: self.view.set_option(attr, v))
            tb.addWidget(cb)
            return cb

        self.cb_h = check("Hydrogens", "show_h", True)
        self.cb_labels = check("Labels", "show_labels", False)
        self.cb_cell = check("Unit cell", "show_cell", False)
        self.cb_poly = check("Polyhedron", "show_polyhedron", True)
        self.cb_dark = check("Dark background", "dark_bg", False)
        tb.addSeparator()

        tb.addWidget(QLabel(" View along: "))
        for ax in "abc":
            b = QPushButton(ax)
            b.setFixedWidth(28)
            b.setToolTip(f"View along the {ax} axis")
            b.clicked.connect(lambda _=False, a=ax: self.view.view_along(a))
            tb.addWidget(b)
        reset = QPushButton("Reset view")
        reset.clicked.connect(self.view.reset_view)
        tb.addWidget(reset)

    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        m.addAction(self.act_open)
        m.addSeparator()
        quit_act = QAction("Exit", self)
        quit_act.setShortcut(QKeySequence.Quit)
        quit_act.triggered.connect(self.close)
        m.addAction(quit_act)
        # a clickable Compare entry in the menu bar (between File and View)
        self.act_compare = QAction("Compare", self)
        self.act_compare.setToolTip("Start a new comparison. While the Compare window is open, every complex "
                                    "you click is added to it.")
        self.act_compare.triggered.connect(self._open_compare)
        self.menuBar().addAction(self.act_compare)
        self.menu_view = self.menuBar().addMenu("&View")
        h = self.menuBar().addMenu("&Help")
        controls = QAction("Controls", self)
        controls.triggered.connect(self._show_help)
        h.addAction(controls)
        about = QAction(f"About {__app_name__}", self)
        about.triggered.connect(self._show_about)
        h.addAction(about)

    def _build_side_panel(self):
        self.side = QWidget()
        lay = QVBoxLayout(self.side)
        lay.setContentsMargins(6, 6, 6, 6)

        lay.addWidget(QLabel("<b>Dataset information</b>"))
        self.info = QTextBrowser()
        self.info.setMaximumHeight(230)
        lay.addWidget(self.info)

        lay.addWidget(QLabel("<b>Inner coordination sphere</b>"))
        self.sphere_label = QLabel("Click on a complex in the view to select its metal centre")
        self.sphere_label.setWordWrap(True)
        lay.addWidget(self.sphere_label)
        self.tabs = QTabWidget()
        lay.addWidget(self.tabs, 1)
        lay_sphere = QVBoxLayout()

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Atom", "x (Å)", "y (Å)", "z (Å)", "d(M–L) Å", "Symmetry"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        lay_sphere.addWidget(self.table, 1)

        row = QHBoxLayout()
        self.btn_copy = QPushButton("Copy coordinates")
        self.btn_copy.setToolTip("Copy as tab-separated text (paste directly into Excel)")
        self.btn_copy.clicked.connect(self._copy_sphere)
        self.btn_clear = QPushButton("Clear selection")
        self.btn_clear.clicked.connect(lambda: self._on_atom_clicked(-1))
        row.addWidget(self.btn_copy)
        row.addWidget(self.btn_clear)
        lay_sphere.addLayout(row)

        tab_sphere = QWidget()
        tab_sphere.setLayout(lay_sphere)
        self.tabs.addTab(tab_sphere, "Sphere")

        # ---- point-group symmetry ----
        tab_sym = QWidget()
        ls = QVBoxLayout(tab_sym)
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Tolerance (Å):"))
        self.edit_tol = QLineEdit(f"{sym.DEFAULT_TOL:g}")
        self.edit_tol.setValidator(QDoubleValidator(1e-9, 1.0, 9))
        self.edit_tol.setMaximumWidth(90)
        self.edit_tol.setToolTip("Maximum atomic displacement allowed when testing a symmetry operation;\n"
                                 "also used as the threshold for Q = 0 and for Q relations.")
        self.edit_tol.editingFinished.connect(self._update_symmetry)
        row2.addWidget(self.edit_tol)
        self.cb_symel = QCheckBox("Show elements")
        self.cb_symel.setChecked(True)
        self.cb_symel.toggled.connect(lambda v: self.view.set_option("show_symmetry", v))
        row2.addWidget(self.cb_symel)
        row2.addStretch(1)
        ls.addLayout(row2)
        self.sym_label = QLabel("—")
        self.sym_label.setTextFormat(Qt.RichText)
        self.sym_label.setWordWrap(True)
        ls.addWidget(self.sym_label)
        self.sym_text = QTextBrowser()
        ls.addWidget(self.sym_text, 1)
        self.tabs.addTab(tab_sym, "Symmetry")

        # ---- van Vleck ----
        tab_vv = QWidget()
        lv = QVBoxLayout(tab_vv)
        row3 = QHBoxLayout()
        self.cb_frame = QCheckBox("Show frame axes")
        self.cb_frame.setChecked(True)
        self.cb_frame.toggled.connect(lambda v: self.view.set_option("show_frame", v))
        row3.addWidget(self.cb_frame)
        row3.addStretch(1)
        self.btn_copy_q = QPushButton("Copy Q")
        self.btn_copy_q.setToolTip("Copy Q1–Q15, constraints and standard-frame coordinates as tab-separated text")
        self.btn_copy_q.clicked.connect(self._copy_q)
        row3.addWidget(self.btn_copy_q)
        self.btn_export_q = QPushButton("Export .txt")
        self.btn_export_q.setToolTip("Save Q1–Q15 of this molecule as a tab-separated .txt file")
        self.btn_export_q.clicked.connect(self._export_q)
        row3.addWidget(self.btn_export_q)
        lv.addLayout(row3)
        self.mass_edit = MassEditor()
        self.mass_edit.changed.connect(self.set_masses)
        lv.addWidget(self.mass_edit)
        self.vv_label = QLabel("—")
        self.vv_label.setTextFormat(Qt.RichText)
        self.vv_label.setWordWrap(True)
        lv.addWidget(self.vv_label)
        self.vv_table = QTableWidget(15, 5)
        self.vv_table.setHorizontalHeaderLabels(["Mode", "Irrep", "Q (Å)", "Symmetry constraint", "Check"])
        self.vv_table.verticalHeader().setVisible(False)
        self.vv_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.vv_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.vv_table.horizontalHeader().setStretchLastSection(True)
        lv.addWidget(self.vv_table, 3)
        self.vv_text = QTextBrowser()
        lv.addWidget(self.vv_text, 2)
        self.tabs.addTab(tab_vv, "van Vleck")
        self._set_sphere_buttons(False)

    def _build_navigator(self):
        """Workspace navigator: files at the top level, datasets below; files can be removed."""
        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(2, 2, 2, 2)
        self.nav = QTreeWidget()
        self.nav.setColumnCount(5)
        self.nav.setHeaderLabels(["Data block", "Compound", "Space group", "T (K)", "CCDC"])
        self.nav.setSelectionMode(QAbstractItemView.SingleSelection)
        self.nav.setRootIsDecorated(True)
        self.nav.setUniformRowHeights(True)
        self.nav.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.nav.currentItemChanged.connect(self._on_nav_item_changed)
        self.nav.setContextMenuPolicy(Qt.CustomContextMenu)
        self.nav.customContextMenuRequested.connect(self._nav_context_menu)
        lay.addWidget(self.nav, 1)
        row = QHBoxLayout()
        b_open = QPushButton("Add files…")
        b_open.clicked.connect(self.open_dialog)
        self.btn_remove_file = QPushButton("Remove file")
        self.btn_remove_file.setToolTip("Remove the selected file from the workspace (Delete key)")
        self.btn_remove_file.clicked.connect(self._remove_selected_file)
        b_clear = QPushButton("Clear")
        b_clear.setToolTip("Remove all files from the workspace")
        b_clear.clicked.connect(self._clear_workspace)
        for b_ in (b_open, self.btn_remove_file, b_clear):
            row.addWidget(b_)
        lay.addLayout(row)
        del_act = QAction(self.nav)
        del_act.setShortcut(QKeySequence.Delete)
        del_act.setShortcutContext(Qt.WidgetShortcut)
        del_act.triggered.connect(self._remove_selected_file)
        self.nav.addAction(del_act)

        dock = QDockWidget("Workspace", self)
        dock.setObjectName("navigator")
        dock.setWidget(box)
        dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea | Qt.BottomDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.resizeDocks([dock], [360], Qt.Horizontal)
        self.menu_view.addAction(dock.toggleViewAction())

    def _rebuild_workspace(self):
        """Rebuild the flat dataset list, the navigator tree and the dataset combo box from self.files."""
        self.structures = [st for f in self.files for st in f["structures"]]
        multi = len(self.files) > 1
        self.nav.blockSignals(True)
        self.nav.clear()
        bold = self.nav.font()
        bold.setBold(True)
        idx = 0
        for f in self.files:
            top = QTreeWidgetItem([os.path.basename(f["path"]), f"{len(f['structures'])} dataset(s)"])
            top.setToolTip(0, f["path"])
            top.setFont(0, bold)
            top.setData(0, Qt.UserRole, ("file", f["path"]))
            self.nav.addTopLevelItem(top)
            for st in f["structures"]:
                vals = [st.short_name, st.short_compound, st.space_group,
                        "" if st.temperature is None else f"{st.temperature:g}", st.ccdc]
                it = QTreeWidgetItem(vals)
                it.setData(0, Qt.UserRole, ("ds", idx))
                for c in range(5):
                    it.setToolTip(c, f"{st.name}\n{st.compound_name}\n{st.formula}")
                top.addChild(it)
                idx += 1
            top.setExpanded(True)
        self.nav.blockSignals(False)

        self.combo_data.blockSignals(True)
        self.combo_data.clear()
        k = 0
        for f in self.files:
            prefix = f"{os.path.basename(f['path'])} › " if multi else ""
            for j, st in enumerate(f["structures"]):
                self.combo_data.addItem(f"{prefix}{j + 1}. {st.display_name()}")
                self.combo_data.setItemData(k, f"{f['path']}\n{st.compound_name or st.name}", Qt.ToolTipRole)
                k += 1
        self.combo_data.blockSignals(False)
        n = len(self.files)
        self.setWindowTitle("vanVleck" if n == 0 else
                            f"vanVleck — {os.path.basename(self.files[0]['path'])}" if n == 1 else
                            f"vanVleck — {n} files")
        self.btn_remove_file.setEnabled(n > 0)

    def _nav_item_for_index(self, idx: int):
        for i in range(self.nav.topLevelItemCount()):
            top = self.nav.topLevelItem(i)
            for j in range(top.childCount()):
                if top.child(j).data(0, Qt.UserRole) == ("ds", idx):
                    return top.child(j)
        return None

    def _on_nav_item_changed(self, item, _prev=None):
        if item is None:
            return
        kind, val = item.data(0, Qt.UserRole)
        if kind == "ds":
            self._on_dataset_changed(val)

    def _file_of_item(self, item) -> str | None:
        if item is None:
            return None
        kind, val = item.data(0, Qt.UserRole)
        if kind == "file":
            return val
        return item.parent().data(0, Qt.UserRole)[1] if item.parent() else None

    def _nav_context_menu(self, pos):
        item = self.nav.itemAt(pos)
        path = self._file_of_item(item)
        if path is None:
            return
        menu = QMenu(self)
        act = menu.addAction(f"Remove “{os.path.basename(path)}” from workspace")
        if menu.exec(self.nav.viewport().mapToGlobal(pos)) == act:
            self.remove_file(path)

    def _remove_selected_file(self):
        path = self._file_of_item(self.nav.currentItem())
        if path is None and len(self.files) == 1:
            path = self.files[0]["path"]
        if path is not None:
            self.remove_file(path)

    def remove_file(self, path: str):
        cur = self.current
        self.files = [f for f in self.files if f["path"] != path]
        self._rebuild_workspace()
        if cur is not None and cur in self.structures:
            self._on_dataset_changed(self.structures.index(cur))
        elif self.structures:
            self._on_dataset_changed(0)
        else:
            self._clear_view()
        self.statusBar().showMessage(f"Removed {os.path.basename(path)} from the workspace", 4000)

    def _clear_workspace(self):
        if not self.files:
            return
        self.files = []
        self._rebuild_workspace()
        self._clear_view()
        self.statusBar().showMessage("Workspace cleared", 4000)

    def _clear_view(self):
        self.current = None
        self.view.set_structure(None)
        self.info.clear()
        self._show_sphere(None)

    # ------------------------------------------------------------ file loading
    def open_dialog(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Add structure files", "", FILE_FILTER)
        self.load_many(paths)

    def load_many(self, paths):
        first = None
        for p in paths:
            idx = self.load(p, select=False)
            if idx is not None and first is None:
                first = idx
        if first is not None:
            self._on_dataset_changed(first)

    def load(self, path: str, select: bool = True):
        """Add a file to the workspace (re-reads it if already present); return the index of its first dataset."""
        path = os.path.abspath(path)
        try:
            QGuiApplication.setOverrideCursor(Qt.WaitCursor)
            structures = load_file(path)
        except Exception as exc:  # noqa: BLE001
            QGuiApplication.restoreOverrideCursor()
            QMessageBox.critical(self, "Read error", f"Could not read file:\n{path}\n\n{exc}")
            return None
        QGuiApplication.restoreOverrideCursor()
        if not structures:
            QMessageBox.warning(self, "No structure data",
                                f"No data block with atomic coordinates was found in\n{path}")
            return None
        key = os.path.normcase(path)
        existing = [i for i, f in enumerate(self.files) if os.path.normcase(f["path"]) == key]
        if existing:
            self.files[existing[0]] = {"path": path, "structures": structures}
            msg = f"Reloaded {os.path.basename(path)}: {len(structures)} dataset(s)"
        else:
            self.files.append({"path": path, "structures": structures})
            msg = f"Added {os.path.basename(path)}: {len(structures)} dataset(s) — {len(self.files)} file(s) in workspace"
        cur = self.current
        self._rebuild_workspace()
        first = self.structures.index(structures[0])
        if select:
            self._on_dataset_changed(first)
        elif cur is not None and cur in self.structures:
            self._on_dataset_changed(self.structures.index(cur))
        self.statusBar().showMessage(msg, 5000)
        return first

    def _on_dataset_changed(self, idx: int):
        if not (0 <= idx < len(self.structures)):
            return
        self.current = self.structures[idx]
        if self.combo_data.currentIndex() != idx:
            self.combo_data.blockSignals(True)
            self.combo_data.setCurrentIndex(idx)
            self.combo_data.blockSignals(False)
        it = self._nav_item_for_index(idx)
        if it is not None and self.nav.currentItem() is not it:
            self.nav.blockSignals(True)
            self.nav.setCurrentItem(it)
            self.nav.blockSignals(False)
        self.view.set_structure(self.current)
        self._update_info()
        self._show_sphere(None)

    def _update_info(self):
        st = self.current
        rows = [("Data block", st.name), ("Source", f"{os.path.basename(st.source_file)} ({st.source_format})")]
        if st.compound_name:
            rows.append(("Compound", st.compound_name))
        if st.formula:
            rows.append(("Formula", st.formula))
        if st.space_group:
            rows.append(("Space group", st.space_group))
        if st.cell:
            c = st.cell
            rows.append(("Cell", f"a={c.a:.4f}, b={c.b:.4f}, c={c.c:.4f} Å<br>"
                                 f"α={c.alpha:.3f}°, β={c.beta:.3f}°, γ={c.gamma:.3f}°<br>"
                                 f"V={c.volume:.1f} Å³"))
        if st.temperature is not None:
            rows.append(("Temperature", f"{st.temperature:g} K"))
        metals = st.metal_indices()
        rows.append(("Atoms shown", f"{len(st.atoms)} ({st.n_molecules} molecules/ions)"))
        rows.append(("Metal centres", ", ".join(st.atoms[i].label for i in metals) or "none"))
        html = "<table cellspacing='2'>" + "".join(
            f"<tr><td valign='top'><b>{k}</b>&nbsp;</td><td>{v}</td></tr>" for k, v in rows) + "</table>"
        for n in st.notes:
            html += f"<p style='color:#b05000'>⚠ {n}</p>"
        self.info.setHtml(html)

    # ------------------------------------------------------------ interaction
    def _on_atom_hovered(self, idx: int):
        if idx < 0 or self.current is None:
            self.statusBar().clearMessage()
            return
        a = self.current.atoms[idx]
        msg = f"{a.label} ({a.element})   x={a.cart[0]:.4f}  y={a.cart[1]:.4f}  z={a.cart[2]:.4f} Å"
        if a.frac is not None:
            msg += f"   frac ({a.frac[0]:.4f}, {a.frac[1]:.4f}, {a.frac[2]:.4f})   [{a.symop}]"
        if a.occupancy < 1:
            msg += f"   occupancy {a.occupancy:g}"
        self.statusBar().showMessage(msg)

    def _on_atom_clicked(self, idx: int):
        if self.current is None:
            return
        if idx < 0:
            self._show_sphere(None)
            return
        sphere = sphere_for_clicked_atom(self.current, idx)
        if sphere is None:
            self.statusBar().showMessage(f"The molecule/ion containing {self.current.atoms[idx].label} has no metal centre", 4000)
            self._show_sphere(None)
            return
        self._show_sphere(sphere)

    def _show_sphere(self, sphere: CoordinationSphere | None):
        self.current_sphere = sphere
        self.table.setRowCount(0)
        if sphere is None:
            self.view.clear_highlight()
            self.sphere_label.setText("Click on a complex in the view to select its metal centre")
            self._set_sphere_buttons(False)
            self._update_symmetry()
            return
        self.view.set_highlight(sphere.metal_index, sphere.ligand_indices, sphere.cis_faces())
        m = sphere.metal
        donors = "".join(sorted({a.element for a in sphere.ligands}))
        dp = dist.compute(sphere.coordinates()) if len(sphere.ligands) == 6 else None
        text = (f"Metal centre <b>{m.label}</b>, {len(sphere.ligands)} donor atoms; "
                f"mean d({m.element}–{donors}) = <b>{sphere.mean_distance:.4f} Å</b>")
        if dp is not None:
            text += (f"<br>ζ = <b>{dp.zeta:.4f} Å</b> &nbsp; Δ = <b>{dist.sci(dp.delta)}</b> &nbsp; "
                     f"Σ = <b>{dp.sigma:.1f}°</b> &nbsp; Θ = <b>{dp.theta:.1f}°</b>"
                     f"<br>Θ<sub>min</sub> = {dp.theta_min:.1f}° &nbsp; Θ<sub>max</sub> = {dp.theta_max:.1f}° &nbsp; "
                     f"V<sub>oct</sub> = {dp.volume:.3f} Å³")
        if sphere.suspicious:
            text += "<br><span style='color:#c00'>⚠ Fewer than 6 donors or an M–L distance > 3 Å — please check</span>"
        self.sphere_label.setText(text)

        entries = [(sphere.metal_index, None)] + list(zip(sphere.ligand_indices, sphere.distances))
        self.table.setRowCount(len(entries))
        for r, (idx, d) in enumerate(entries):
            a = self.current.atoms[idx]
            x, y, z = (float(v) + 0.0 for v in a.cart.round(4))
            vals = [a.label, f"{x:.4f}", f"{y:.4f}", f"{z:.4f}",
                    "—" if d is None else f"{d:.4f}", a.symop]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if 1 <= c <= 4:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if a.frac is not None:
                    it.setToolTip(f"Fractional: {a.frac[0]:.5f}, {a.frac[1]:.5f}, {a.frac[2]:.5f}")
                if d is None:
                    f = it.font()
                    f.setBold(True)
                    it.setFont(f)
                self.table.setItem(r, c, it)
        self._set_sphere_buttons(True)
        self._update_symmetry()
        self.statusBar().showMessage(f"Selected inner coordination sphere of {m.label}", 4000)

    def _tolerance(self) -> float:
        try:
            t = float(self.edit_tol.text())
            return t if t > 0 else sym.DEFAULT_TOL
        except ValueError:
            return sym.DEFAULT_TOL

    def _update_symmetry(self):
        sphere = self.current_sphere
        if sphere is None:
            self.current_symmetry = None
            self.sym_label.setText("—")
            self.sym_text.clear()
            self.view.set_symmetry_elements([])
            self._update_vanvleck()
            return
        QGuiApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            elems = [sphere.metal.element] + [a.element for a in sphere.ligands]
            X = sphere.coordinates()
            tol = self._tolerance()
            res = sym.analyze(X, elems, tol)
            scan = sym.symmetry_vs_tolerance(X, elems)
        finally:
            QGuiApplication.restoreOverrideCursor()
        self.current_symmetry = res
        ok = ("<span style='color:#080'>✔ standard ML<sub>6</sub> subgroup</span>" if res.standard else
              "<span style='color:#c00'>✘ outside the 19 standard subgroups</span>")
        self.sym_label.setText(f"<span style='font-size:16pt'><b>{sym.html(res.label)}</b></span>"
                               f"&nbsp;&nbsp;order {res.order} &nbsp; {ok}<br>tolerance {tol:g} Å, "
                               f"max deviation {res.max_deviation:.2e} Å")
        lines = ["<b>Symmetry elements</b> (Cartesian frame of the structure, origin at metal):<br>"]
        lines += [e + "<br>" for e in res.elements_summary()]
        for n in res.notes:
            lines.append(f"<p style='color:#b05000'>⚠ {n}</p>")
        lines.append("<b>Pseudo-symmetry vs tolerance</b><br>")
        lines.append(" &nbsp;·&nbsp; ".join(f"{t:g} Å: {sym.html(lab)}" for t, lab in scan))
        self.sym_text.setHtml("".join(lines))
        m = sphere.metal.cart
        self.view.set_symmetry_elements(
            [(o.kind, o.order, o.axis, m) for o in res.ops if o.kind in ("C", "σ", "S")])
        self._update_vanvleck()

    def _update_vanvleck(self):
        sphere, res = self.current_sphere, self.current_symmetry
        self.current_vanvleck = None
        self.view.set_frame(None, None)
        self.vv_table.clearContents()
        self.btn_copy_q.setEnabled(False)
        self.btn_export_q.setEnabled(False)
        if sphere is None or res is None:
            self.vv_label.setText("—")
            self.vv_text.clear()
            return
        if res.label not in vv.STANDARD_SUBGROUPS:
            self.vv_label.setText(f"No standard frame defined for {sym.html(res.label)}.")
            self.vv_text.clear()
            return
        elems = [sphere.metal.element] + [a.element for a in sphere.ligands]
        X = sphere.coordinates()
        tol = self._tolerance()
        QGuiApplication.setOverrideCursor(Qt.WaitCursor)
        # while Compare is open, the molecule is added there and shown with the frame used in Compare
        cmp_mode = self.compare_dialog is not None and self.compare_dialog.isVisible()
        try:
            r = None
            if cmp_mode:
                entry, added = self.compare_dialog.add(sphere, res, self.current.short_name,
                                                       self.current.temperature, tol)
                self._update_compare_button()
                r = entry.result
            if r is None:
                m0, m = self.mass_edit.values()
                r = vv.compute(X, elems, res, zero_tol=tol, m_metal=m0, m_lig=m)
        except Exception as exc:  # noqa: BLE001
            self.vv_label.setText(f"<span style='color:#c00'>Frame selection failed: {exc}</span>")
            self.vv_text.clear()
            return
        finally:
            QGuiApplication.restoreOverrideCursor()
        self.current_vanvleck = r
        n_bad = r.checks.count("violated")
        status = ("<span style='color:#080'>✔ all zero / non-zero relations satisfied</span>" if n_bad == 0 else
                  f"<span style='color:#c00'>✘ {n_bad} relation(s) violated</span>")
        cmp_txt = "<br>Frame and Q as used in the Compare window." if cmp_mode else ""
        self.vv_label.setText(f"Standard frame for <b>{sym.html(r.label)}</b> &nbsp; "
                              f"|Q| (Q<sub>2</sub>–Q<sub>15</sub>) = <b>{r.norm:.4f} Å</b> &nbsp; {status}{cmp_txt}")
        green, red, amber = QColor(0, 130, 0), QColor(200, 0, 0), QColor(190, 110, 0)
        for k in range(15):
            c = r.constraints[k]
            chk = r.checks[k]
            sym_txt = {"ok": "✔", "violated": "✘", "accidental-zero": "0 (allowed ≠ 0)"}[chk]
            qv = 0.0 if abs(r.Q[k]) < 5e-7 else r.Q[k]
            vals = [f"Q{k + 1}", vv.IRREPS[k], f"{qv:+.6f}", c.text(), sym_txt]
            for col, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if col == 2:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    if abs(r.Q[k]) < tol:
                        it.setForeground(QColor(150, 150, 150))
                if col == 4:
                    it.setForeground(green if chk == "ok" else red if chk == "violated" else amber)
                if col == 3 and c.zero:
                    it.setForeground(QColor(120, 120, 120))
                self.vv_table.setItem(k, col, it)
        names = ["+x", "+y", "+z", "−x", "−y", "−z"]
        lig = sphere.ligands
        html = ["<b>Ligand positions in the standard frame</b> (origin at metal):<br>"
                "<table cellspacing='3'><tr><th>site</th><th>atom</th><th>x</th><th>y</th><th>z</th></tr>"]
        for k in range(6):
            a = lig[r.ligand_at_site[k]]
            x, y, z = (0.0 if abs(v) < 5e-7 else v for v in r.coords_std[k + 1])
            html.append(f"<tr><td>{k + 1} ({names[k]})</td><td>{a.label}</td><td>{x:+.5f}</td>"
                        f"<td>{y:+.5f}</td><td>{z:+.5f}</td></tr>")
        html.append("</table>")
        html.append(f"<p>Alignment misfit Σ|R·u<sub>i</sub> − e<sub>i</sub>|² = {r.misfit:.5f}<br>"
                    f"Masses: metal M₀ = {r.masses[0]:g}, ligand atoms M = {r.masses[1]:g}</p>")
        Rm = r.R
        html.append("<p><b>Rotation matrix</b> (structure → standard frame):<br><tt>" +
                    "<br>".join("&nbsp;".join(f"{v:+.6f}" for v in row) for row in Rm) + "</tt></p>")
        for n in r.notes:
            html.append(f"<p style='color:#b05000'>⚠ {n}</p>")
        self.vv_text.setHtml("".join(html))
        self.btn_copy_q.setEnabled(True)
        self.btn_export_q.setEnabled(True)
        self.view.set_frame(Rm, sphere.metal.cart)
        if cmp_mode:                    # no separate bar chart while comparing
            return
        # bar-chart pop-up
        if self.qdialog is None:
            self.qdialog = QChartDialog(self)
        T = f", {self.current.temperature:g} K" if self.current.temperature is not None else ""
        title = (f"{self.current.short_name}{T} — {sphere.metal.label}, {sym.pretty(r.label)}"
                 f"     |Q| = {r.norm:.4f} Å")
        n_bad = r.checks.count("violated")
        info = (f"Zero threshold {tol:g} Å.  M₀ = {r.masses[0]:g}, M = {r.masses[1]:g}.  " +
                ("All symmetry relations satisfied." if n_bad == 0 else f"{n_bad} symmetry relation(s) violated."))
        dp = dist.compute(X)
        sub = (f"⟨d⟩ = {dp.d_mean:.4f} Å     ζ = {dp.zeta:.4f} Å     Δ = {dist.sci(dp.delta)}     "
               f"Σ = {dp.sigma:.1f}°     Θ = {dp.theta:.1f}°\n"
               f"Θmin = {dp.theta_min:.1f}°     Θmax = {dp.theta_max:.1f}°     V(oct) = {dp.volume:.3f} Å³")
        self.qdialog.show_result(r.Q, vv.IRREPS, tol, title, info, name=self._q_name(), subtitle=sub)

    # ------------------------------------------------------------ compare
    def _ensure_compare(self) -> CompareDialog:
        if self.compare_dialog is None:
            self.compare_dialog = CompareDialog(self)
            self.compare_dialog.changed.connect(self._update_compare_button)
            self.compare_dialog.set_masses(*self.mass_edit.values())
            self.compare_dialog.mass_edit.changed.connect(self.set_masses)
        return self.compare_dialog

    def _update_compare_button(self):
        n = len(self.compare_dialog.entries) if self.compare_dialog else 0
        self.act_compare.setToolTip(f"Open the Compare window ({n} molecules)")

    def _open_compare(self):
        """Each click on Compare starts a new comparison: the previous list is cleared and the
        complex selected in the main window, if any, becomes the first entry."""
        d = self._ensure_compare()
        d.entries.clear()
        d.set_masses(*self.mass_edit.values())
        d.show()
        d.raise_()
        if self.current_sphere is not None and self.current_symmetry is not None:
            d.add(self.current_sphere, self.current_symmetry, self.current.short_name,
                  self.current.temperature, self._tolerance())
        d.recompute()
        self._update_compare_button()
        if self.current_sphere is not None:     # show the Compare frame in the main window
            self._update_vanvleck()

    def set_masses(self, m0: float, m: float):
        """After the masses change in either window: sync both inputs and recompute all Q."""
        self.mass_edit.set_values(m0, m)
        if self.compare_dialog is not None:
            self.compare_dialog.set_masses(m0, m)
        if self.current_sphere is not None:
            self._update_vanvleck()
        self.statusBar().showMessage(f"Masses set to M₀ = {m0:g}, M = {m:g}; Q recalculated", 4000)

    def _add_to_compare(self):
        """Add the selected complex to the comparison (duplicates are ignored)."""
        if self.current_sphere is None or self.current_symmetry is None:
            return
        d = self._ensure_compare()
        added = d.add(self.current_sphere, self.current_symmetry, self.current.short_name,
                      self.current.temperature, self._tolerance())
        self._update_compare_button()
        if added:
            self.statusBar().showMessage(f"Selected inner coordination sphere of {self.current_sphere.metal.label}"
                                         f" — added to Compare ({len(d.entries)})", 4000)

    def _q_name(self) -> str:
        st, sp = self.current, self.current_sphere
        T = f" {st.temperature:g} K" if st.temperature is not None else ""
        return f"{st.short_name}{T} {sp.metal.label}"

    def _export_q(self):
        r = self.current_vanvleck
        if r is None:
            return
        name = self._q_name()
        p = export_q_dialog(self, [name], [r.Q], f"Q_{name}.txt")
        if p:
            self.statusBar().showMessage(f"Saved {p}", 5000)

    def _copy_q(self):
        r = self.current_vanvleck
        if r is None:
            return
        txt = r.to_text() + "\n\nSite\tAtom\tx\ty\tz\n"
        lig = self.current_sphere.ligands
        txt += f"0\t{self.current_sphere.metal.label}\t0\t0\t0\n"
        for k in range(6):
            x, y, z = r.coords_std[k + 1]
            txt += f"{k + 1}\t{lig[r.ligand_at_site[k]].label}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
        QGuiApplication.clipboard().setText(txt)
        self.statusBar().showMessage("van Vleck coordinates copied to clipboard", 3000)

    def _set_sphere_buttons(self, on: bool):
        self.btn_copy.setEnabled(on)
        self.btn_clear.setEnabled(on)

    def _copy_sphere(self):
        if self.current_sphere:
            QGuiApplication.clipboard().setText(self.current_sphere.to_text())
            self.statusBar().showMessage("Coordinates copied to clipboard", 3000)

    def _show_about(self):
        QMessageBox.about(self, f"About {__app_name__}",
            f"<h3>{__app_name__} {__version__}</h3>"
            "<p>Point-group symmetry and van Vleck coordinates Q<sub>1</sub>–Q<sub>15</sub> of "
            "pseudo-octahedral ML<sub>6</sub> coordination spheres, with OctaDist-type distortion "
            "parameters (ζ, Δ, Σ, Θ).</p>"
            f"<p><a href='{__url__}'>{__url__}</a></p>"
            "<p>Released under the MIT License. If you use this program in published work, "
            "please cite it as described in the README / CITATION.cff.</p>")

    def _show_help(self):
        QMessageBox.information(self, "Controls",
            "• Left-drag: rotate\n• Mouse wheel / right-drag up-down: zoom\n• Middle-drag or Shift+left-drag: pan\n"
            "• Click a complex: highlight the metal centre and its 6 donor atoms\n• Click empty space: clear selection\n"
            "• Drag and drop one or more .cif / .xyz files onto the window; they are added to the workspace\n"
            "• Remove a file from the workspace: select it and press Delete, or right-click it\n"
            "• If a CIF contains several data blocks, switch with the 'Dataset' box or the Workspace tree")

    # ------------------------------------------------------------ drag & drop
    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e):
        paths = [u.toLocalFile() for u in e.mimeData().urls()]
        self.load_many([p for p in paths if p.lower().endswith((".cif", ".xyz"))])
