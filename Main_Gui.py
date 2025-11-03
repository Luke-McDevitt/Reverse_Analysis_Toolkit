# Main_Gui.py
# Tkinter GUI shell for the Spherical Equations Toolbox
# Tab 1 = Mobility calculator (spatial or spherical), scrollable, auto joint rows.
# No external deps beyond stdlib. Tested on Python 3.11+.

import math
import tkinter as tk
from tkinter import ttk, messagebox

# Ensure the folder containing Main_Gui.py is on sys.path
import os, sys
CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from Mobility.MobilityTab import MobilityTab

APP_TITLE = "Spherical Equations Toolbox — Main GUI"
PAD = 8

class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1000x720")

        # Styles
        style = ttk.Style(self)
        style.configure("Headline.TLabel", font=("Segoe UI", 14, "bold"))
        style.configure("Header.TLabel", font=("Segoe UI", 10, "bold"))

        # Notebook with tabs (Mobility first)
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)

        self.tab_mobility = MobilityTab(nb)
        nb.add(self.tab_mobility, text="Mobility")

        # Placeholder tabs you can fill later
        placeholder = ttk.Frame(nb)
        ttk.Label(placeholder, text="(Next: Toolbox tabs go here)").pack(padx=PAD, pady=PAD, anchor="w")
        nb.add(placeholder, text="(More Tabs)")

        # Menubar (simple)
        self._build_menu()

    def _build_menu(self):
        m = tk.Menu(self)
        filem = tk.Menu(m, tearoff=False)
        filem.add_command(label="Exit", command=self.destroy)
        m.add_cascade(label="File", menu=filem)

        helpm = tk.Menu(m, tearoff=False)
        helpm.add_command(label="About", command=lambda: messagebox.showinfo("About", APP_TITLE))
        m.add_cascade(label="Help", menu=helpm)
        self.config(menu=m)


if __name__ == "__main__":
    MainApp().mainloop()
