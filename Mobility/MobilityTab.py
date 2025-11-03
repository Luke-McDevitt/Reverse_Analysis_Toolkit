# Mobility/MobilityTab.py
# Standalone Mobility tab with Save/Load + Examples
import os, csv, glob, math, datetime
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from functions_and_helpers.mobility import compute_mobility, JOINT_FI_DEFAULTS

PAD = 8

def _project_root():
    # default file dialogs start in the folder that contains Main_Gui.py (one level up)
    return os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

def _tab_dir():
    return os.path.abspath(os.path.dirname(__file__))

class ScrollFrame(ttk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.vbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)

        self.inner.bind("<Configure>", self._on_cfg)
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.configure(yscrollcommand=self.vbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.vbar.pack(side="right", fill="y")

        # Mouse wheel (Win/Mac/Linux)
        self.canvas.bind_all("<MouseWheel>", self._on_wheel)
        self.canvas.bind_all("<Button-4>", self._on_wheel)  # Linux up
        self.canvas.bind_all("<Button-5>", self._on_wheel)  # Linux down

    def _on_cfg(self, _e=None):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        self.canvas.itemconfigure(self._win, width=self.canvas.winfo_width())

    def _on_wheel(self, e):
        if getattr(e, "num", None) == 4:
            delta = -120
        elif getattr(e, "num", None) == 5:
            delta = 120
        else:
            delta = -1 * (e.delta if getattr(e, "delta", 0) else 0)
        self.canvas.yview_scroll(int(delta/120), "units")

class MobilityTab(ttk.Frame):
    """Mobility calculator tab with Save/Load and Examples list."""
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # Styles (safe to call repeatedly)
        style = ttk.Style(self)
        style.configure("Headline.TLabel", font=("Segoe UI", 14, "bold"))
        style.configure("Header.TLabel",   font=("Segoe UI", 10, "bold"))

        self.scroll = ScrollFrame(self)
        self.scroll.pack(fill="both", expand=True)
        root = self.scroll.inner

        # Title
        ttk.Label(root, text="Mobility Calculator", style="Headline.TLabel").grid(
            row=0, column=0, columnspan=8, sticky="w", padx=PAD, pady=(PAD, 0)
        )

        # Save/Load/Examples bar (per Program Plan)
        bar = ttk.Frame(root)
        bar.grid(row=1, column=0, columnspan=8, sticky="w", padx=PAD, pady=(4, PAD))
        ttk.Button(bar, text="Save", command=self._on_save).pack(side="left")
        ttk.Button(bar, text="Load", command=self._on_load).pack(side="left", padx=(6,0))
        ttk.Label(bar, text="Examples:").pack(side="left", padx=(12,4))
        self.examples_cb_var = tk.StringVar(value="")
        self.examples_cb = ttk.Combobox(bar, textvariable=self.examples_cb_var, state="readonly", width=42)
        self.examples_cb.pack(side="left")
        ttk.Button(bar, text="Load Example", command=self._on_load_example).pack(side="left", padx=(6,0))

        # Spherical option
        self.is_spherical_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(root, text="This is a spherical loop", variable=self.is_spherical_var,
                        command=self._recompute).grid(row=2, column=0, columnspan=2, sticky="w", padx=PAD, pady=(0, PAD))

        # Counts row
        frm_counts = ttk.Frame(root)
        frm_counts.grid(row=3, column=0, columnspan=8, sticky="ew", padx=PAD, pady=(0, PAD))
        frm_counts.columnconfigure(7, weight=1)

        ttk.Label(frm_counts, text="Number of bodies (n):").grid(row=0, column=0, sticky="w", padx=(0,4))
        self.n_var = tk.IntVar(value=5)
        ttk.Spinbox(frm_counts, from_=1, to=100, textvariable=self.n_var, width=6,
                    command=self._on_counts_change).grid(row=0, column=1, sticky="w", padx=(0,16))

        ttk.Label(frm_counts, text="Number of joints (j):").grid(row=0, column=2, sticky="w", padx=(0,4))
        self.j_var = tk.IntVar(value=5)
        ttk.Spinbox(frm_counts, from_=0, to=100, textvariable=self.j_var, width=6,
                    command=self._on_counts_change).grid(row=0, column=3, sticky="w", padx=(0,16))

        ttk.Label(frm_counts, text="ΔM (special geometry):").grid(row=0, column=4, sticky="e", padx=(0,4))
        self.deltaM_var = tk.DoubleVar(value=0.0)
        ttk.Entry(frm_counts, textvariable=self.deltaM_var, width=8).grid(row=0, column=5, sticky="w")

        self.chain_info = ttk.Label(root, text="", foreground="#666")
        self.chain_info.grid(row=4, column=0, columnspan=8, sticky="w", padx=PAD, pady=(0, PAD))

        # Joint table
        self.joint_rows_frame = ttk.Frame(root)
        self.joint_rows_frame.grid(row=5, column=0, columnspan=8, sticky="ew", padx=PAD, pady=PAD)
        self.joint_rows_frame.columnconfigure(6, weight=1)
        self._build_joint_table_header()
        self.joint_rows = []
        self._rebuild_joint_rows()

        # Quick inserts
        ex_frame = ttk.Frame(root)
        ex_frame.grid(row=6, column=0, columnspan=8, sticky="w", padx=PAD, pady=(0, PAD))
        ttk.Label(ex_frame, text="Quick inserts:").pack(side="left")
        ttk.Button(ex_frame, text="3R-2C (5-link)", command=self._load_ex_3R2C).pack(side="left", padx=(6,0))
        ttk.Button(ex_frame, text="6R closed", command=self._load_ex_6R).pack(side="left", padx=6)
        ttk.Button(ex_frame, text="Spherical 4R", command=self._load_ex_sph4R).pack(side="left", padx=6)

        # Output options
        opt_frame = ttk.LabelFrame(root, text="Output options")
        opt_frame.grid(row=7, column=0, columnspan=8, sticky="ew", padx=PAD, pady=(0, PAD))
        self.show_steps_var = tk.BooleanVar(value=True)
        self.show_verbose_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt_frame, text="Show step-by-step equations", variable=self.show_steps_var,
                        command=self._recompute).pack(side="left", padx=(4, 12))
        ttk.Checkbutton(opt_frame, text="Explain solution verbally", variable=self.show_verbose_var,
                        command=self._recompute).pack(side="left", padx=(0, 12))

        # Results
        out_frame = ttk.LabelFrame(root, text="Results")
        out_frame.grid(row=8, column=0, columnspan=8, sticky="nsew", padx=PAD, pady=(0, PAD))
        out_frame.columnconfigure(0, weight=1)
        out_frame.rowconfigure(0, weight=1)

        self.output_text = tk.Text(out_frame, width=100, height=16, wrap="word")
        self.output_text.grid(row=0, column=0, sticky="nsew")

        btn_row = ttk.Frame(root)
        btn_row.grid(row=9, column=0, columnspan=8, sticky="e", padx=PAD, pady=(0, PAD))
        ttk.Button(btn_row, text="Compute Mobility", command=self._recompute).pack(side="right")
        ttk.Button(btn_row, text="Copy Output", command=self._copy_output).pack(side="right", padx=(0, 8))

        # Build examples dropdown from files in this folder
        self._refresh_examples_list()

        # First compute
        self._recompute()

    # ---------- UI helpers ----------
    def _build_joint_table_header(self):
        hdrs = ["#", "Type", "fi (DOF)", "Connects body a", "to body b"]
        for c, t in enumerate(hdrs):
            ttk.Label(self.joint_rows_frame, text=t, style="Header.TLabel").grid(row=0, column=c, padx=4, pady=(0,4), sticky="w")

    def _rebuild_joint_rows(self):
        if hasattr(self, "joint_rows"):
            for (widgets, _) in self.joint_rows:
                for w in widgets: w.destroy()
            self.joint_rows.clear()
        j = max(0, int(self.j_var.get()))
        for i in range(j):
            r = i + 1
            idx = ttk.Label(self.joint_rows_frame, text=str(i+1))
            idx.grid(row=r, column=0, sticky="w", padx=4, pady=2)

            typ_var = tk.StringVar(value="R")
            cb = ttk.Combobox(self.joint_rows_frame, textvariable=typ_var,
                              values=("R", "P", "C", "B", "PL", "Custom"), width=8, state="readonly")
            cb.grid(row=r, column=1, sticky="w", padx=4, pady=2)

            fi_var = tk.DoubleVar(value=1.0)
            fi_entry = ttk.Entry(self.joint_rows_frame, textvariable=fi_var, width=8)
            fi_entry.grid(row=r, column=2, sticky="w", padx=4, pady=2)

            def _on_sel(_e=None, tv=typ_var, fv=fi_var, e=fi_entry):
                t = tv.get().upper()
                if t in JOINT_FI_DEFAULTS:
                    fv.set(JOINT_FI_DEFAULTS[t]);
                    e.configure(state="normal")
                else:
                    e.configure(state="normal")

            cb.bind("<<ComboboxSelected>>", _on_sel)

            a_var = tk.IntVar(value=max(1, min(i+1, 99)))
            b_var = tk.IntVar(value=max(1, min(i+2, 99)))
            a_entry = ttk.Entry(self.joint_rows_frame, textvariable=a_var, width=10)
            b_entry = ttk.Entry(self.joint_rows_frame, textvariable=b_var, width=10)
            a_entry.grid(row=r, column=3, sticky="w", padx=4, pady=2)
            b_entry.grid(row=r, column=4, sticky="w", padx=4, pady=2)

            widgets = (idx, cb, fi_entry, a_entry, b_entry)
            state = {"typ": typ_var, "fi": fi_var, "a": a_var, "b": b_var}
            self.joint_rows.append((widgets, state))

        n = int(self.n_var.get())
        if j == n:
            self.chain_info.configure(text="Detected single-chain closed loop (n=j). Special forms will be shown.")
        else:
            self.chain_info.configure(text="")

    def _on_counts_change(self):
        try:
            self.n_var.set(int(self.n_var.get()))
            self.j_var.set(int(self.j_var.get()))
        except Exception:
            pass
        self._rebuild_joint_rows()
        self._recompute()

    # ---------- Quick inserts ----------
    def _load_ex_3R2C(self):
        self.is_spherical_var.set(False)
        self.n_var.set(5); self.j_var.set(5); self._rebuild_joint_rows()
        types = ["R","C","R","C","R"]
        for (_, st), t in zip(self.joint_rows, types):
            st["typ"].set(t); st["fi"].set(2.0 if t=="C" else 1.0)
        self._recompute()

    def _load_ex_6R(self):
        self.is_spherical_var.set(False)
        self.n_var.set(6); self.j_var.set(6); self._rebuild_joint_rows()
        for (_, st) in self.joint_rows:
            st["typ"].set("R"); st["fi"].set(1.0)
        self._recompute()

    def _load_ex_sph4R(self):
        self.is_spherical_var.set(True)
        self.n_var.set(4); self.j_var.set(4); self._rebuild_joint_rows()
        for (_, st) in self.joint_rows:
            st["typ"].set("R"); st["fi"].set(1.0)
        self._recompute()

    # ---------- Compute ----------
    @staticmethod
    def _sum_fi(rows):
        vals = []
        for _, st in rows:
            try: vals.append(float(st["fi"].get()))
            except Exception: vals.append(0.0)
        return sum(vals), vals

    def _recompute(self):
        try:
            n = int(self.n_var.get())
            j = int(self.j_var.get())
            deltaM = float(self.deltaM_var.get())
        except Exception:
            messagebox.showerror("Invalid inputs", "Please enter valid numeric values.")
            return

        spherical = self.is_spherical_var.get()
        show_steps = self.show_steps_var.get()
        show_verbose = self.show_verbose_var.get()

        # Collect fi values from the joint rows (we keep types for UI, but only fi are needed here)
        fi_list = []
        for _, st in self.joint_rows:
            try:
                fi_list.append(float(st["fi"].get()))
            except Exception:
                fi_list.append(0.0)

        # Use shared helper
        results = compute_mobility(n=n, fi_list=fi_list, spherical=spherical, deltaM=deltaM)
        sum_fi = results["sum_fi"]
        M_spatial_general = results["spatial_general"]
        M_spatial_closed = results["spatial_closed"]
        M_spherical_general = results["spherical_general"]
        M_spherical_closed = results["spherical_closed"]
        group_num = results["group"]

        # Pick the "primary" readout to match current toggle & single-chain status
        n_eq_j = (j == n)
        if spherical:
            if n_eq_j:
                M_primary = M_spherical_closed
                primary_name = "Spherical (single-chain closed loop)"
                primary_eq_sym = "M = Σfi − 3"
                primary_eq_num = f"M = {sum_fi:.3g} − 3"
            else:
                M_primary = M_spherical_general
                primary_name = "Spherical (general)"
                primary_eq_sym = "M = 3(n − 1) − Σ(3 − fi)"
                primary_eq_num = f"M = 3({n} − 1) − Σ(3 − fi)"
        else:
            if n_eq_j:
                M_primary = M_spatial_closed
                primary_name = "Spatial (single-chain closed loop)"
                primary_eq_sym = "M = Σfi − 6"
                primary_eq_num = f"M = {sum_fi:.3g} − 6"
            else:
                M_primary = M_spatial_general
                primary_name = "Spatial (general)"
                primary_eq_sym = "M = 6(n − 1) − Σ(6 − fi)"
                primary_eq_num = f"M = 6({n} − 1) − Σ(6 − fi)"

        # ---------- Render ----------
        self.output_text.configure(state="normal")
        self.output_text.delete("1.0", "end")

        self.output_text.insert("end", f"Primary result ({primary_name}):\n  {primary_eq_sym}\n")
        if show_steps:
            if spherical and (not n_eq_j):
                sigma_terms = " + ".join([f"(3 − {fi:g})" for fi in fi_list]) or "0"
                self.output_text.insert("end", f"  {primary_eq_num} = 3*({n - 1}) − ({sigma_terms})\n")
            elif (not spherical) and (not n_eq_j):
                sigma_terms = " + ".join([f"(6 − {fi:g})" for fi in fi_list]) or "0"
                self.output_text.insert("end", f"  {primary_eq_num} = 6*({n - 1}) − ({sigma_terms})\n")
            else:
                self.output_text.insert("end", f"  {primary_eq_num}\n")

        self.output_text.insert("end", f"  ⇒ Mobility M = {M_primary:.6g}\n\n")

        # Show all forms for reference
        self.output_text.insert("end", "Reference (all forms):\n")
        self.output_text.insert("end", f"  Spatial general:        M = 6(n−1) − Σ(6−fi)  ⇒ {M_spatial_general:.6g}\n")
        if n_eq_j:
            self.output_text.insert("end", f"  Spatial closed:         M = Σfi − 6          ⇒ {M_spatial_closed:.6g}\n")
        self.output_text.insert("end", f"  Spherical general:      M = 3(n−1) − Σ(3−fi)  ⇒ {M_spherical_general:.6g}\n")
        if n_eq_j:
            self.output_text.insert("end",
                                    f"  Spherical closed:       M = Σfi − 3          ⇒ {M_spherical_closed:.6g}\n")
        if (group_num is not None) and (not spherical):
            self.output_text.insert("end", f"  Group (equiv spherical):G = Σfi − 3          ⇒ {group_num:.6g}\n")

        if show_verbose:
            self.output_text.insert("end", "\nExplanation:\n")
            if spherical:
                self.output_text.insert("end",
                                        "  Each link has 3 DOF: M = 3(n−1) − Σ(3−fi). For a single-chain loop (n=j), M = Σfi − 3.\n")
            else:
                self.output_text.insert("end",
                                        "  Each link has 6 DOF: M = 6(n−1) − Σ(6−fi). For a single-chain loop (n=j), M = Σfi − 6. Equivalent spherical mobility (Group) is G = Σfi − 3.\n")
            self.output_text.insert("end", "  ΔM accounts for special geometry.\n")

        self.output_text.insert("end", "\nInputs summary:\n")
        self.output_text.insert("end",
                                f"  n = {n}, j = {j}, spherical={bool(self.is_spherical_var.get())}, ΔM = {deltaM}\n")
        self.output_text.insert("end",
                                f"  fi list (per joint): {', '.join(f'{fi:g}' for fi in fi_list) if fi_list else '(none)'}\n")
        self.output_text.configure(state="disabled")

    def _copy_output(self):
        txt = self.output_text.get("1.0","end-1c")
        self.clipboard_clear(); self.clipboard_append(txt)
        messagebox.showinfo("Copied", "Results copied to clipboard.")

    # ---------- Save/Load/Examples ----------
    def _state_to_dict(self):
        d = {
            "spherical": int(self.is_spherical_var.get()),
            "n": int(self.n_var.get()),
            "j": int(self.j_var.get()),
            "deltaM": float(self.deltaM_var.get()),
            "show_steps": int(self.show_steps_var.get()),
            "show_verbose": int(self.show_verbose_var.get()),
        }
        for i, (_, st) in enumerate(self.joint_rows):
            d[f"joint[{i}].type"] = st["typ"].get()
            d[f"joint[{i}].fi"]   = st["fi"].get()
            d[f"joint[{i}].a"]    = st["a"].get()
            d[f"joint[{i}].b"]    = st["b"].get()
        return d

    def _apply_state_from_dict(self, d):
        sph = bool(int(d.get("spherical", 0)))
        n   = int(d.get("n", 1))
        # j from explicit or from highest joint index
        joint_idxs = sorted({int(k.split("[",1)[1].split("]",1)[0]) for k in d if k.startswith("joint[")})
        j = int(d.get("j", len(joint_idxs)))
        if joint_idxs: j = max(j, 1 + max(joint_idxs))

        self.is_spherical_var.set(sph)
        self.n_var.set(n); self.j_var.set(j)
        self._rebuild_joint_rows()

        try: self.deltaM_var.set(float(d.get("deltaM", 0.0)))
        except: pass
        try: self.show_steps_var.set(bool(int(d.get("show_steps", 1))))
        except: pass
        try: self.show_verbose_var.set(bool(int(d.get("show_verbose", 1))))
        except: pass

        for i in range(min(j, len(self.joint_rows))):
            _, st = self.joint_rows[i]
            st["typ"].set(d.get(f"joint[{i}].type", st["typ"].get()))
            try: st["fi"].set(float(d.get(f"joint[{i}].fi", st["fi"].get())))
            except: pass
            try: st["a"].set(int(d.get(f"joint[{i}].a", st["a"].get())))
            except: pass
            try: st["b"].set(int(d.get(f"joint[{i}].b", st["b"].get())))
            except: pass

        self._recompute()

    def _on_save(self):
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        default = f"{now}_Mobility.csv"
        path = filedialog.asksaveasfilename(
            title="Save Mobility state",
            initialdir=_project_root(),
            initialfile=default,
            defaultextension=".csv",
            filetypes=[("CSV files","*.csv")]
        )
        if not path: return
        d = self._state_to_dict()
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f); w.writerow(["key","value"])
                for k,v in d.items(): w.writerow([k, v])
            messagebox.showinfo("Saved", f"Saved state to:\n{path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))

    def _on_load(self):
        path = filedialog.askopenfilename(
            title="Load Mobility state",
            initialdir=_project_root(),
            filetypes=[("CSV files","*.csv")]
        )
        if not path: return
        try:
            with open(path, "r", newline="", encoding="utf-8") as f:
                r = csv.reader(f); rows = list(r)
            if rows and rows[0][:2] == ["key","value"]: rows = rows[1:]
            d = {k:v for (k,v) in rows if k}
            self._apply_state_from_dict(d)
            messagebox.showinfo("Loaded", f"Loaded state from:\n{path}")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))

    def _refresh_examples_list(self):
        pattern = os.path.join(_tab_dir(), "Example_*_Mobility.csv")
        files = sorted(glob.glob(pattern))
        self._example_map = {}
        display_names = []
        for p in files:
            base = os.path.basename(p)
            name = base.removeprefix("Example_").removesuffix("_Mobility.csv")
            display_names.append(name)
            self._example_map[name] = p
        self.examples_cb["values"] = display_names
        if display_names:
            self.examples_cb_var.set(display_names[0])

    def _on_load_example(self):
        name = self.examples_cb_var.get()
        path = self._example_map.get(name, "")
        if not path:
            messagebox.showwarning("No example selected", "Please choose an example first.")
            return
        try:
            with open(path, "r", newline="", encoding="utf-8") as f:
                r = csv.reader(f); rows = list(r)
            if rows and rows[0][:2] == ["key","value"]: rows = rows[1:]
            d = {k:v for (k,v) in rows if k}
            self._apply_state_from_dict(d)
            messagebox.showinfo("Example loaded", f"Loaded example:\n{name}")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))
