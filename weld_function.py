
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import math
import numpy as np


class WeldingSimulation:
    # 初期化関数
    def __init__(self, x3, x4, A1, lambda_1, mu, sigma, lambda_3, C, alpha, amp):
        # 各パラメータをクラス内部変数として保持
        self.A1 = A1
        self.lambda_1 = lambda_1
        self.mu = mu
        self.sigma = sigma
        self.lambda_3 = lambda_3
        self.C = C
        self.alpha = alpha
        self.amp = amp

        # 各変数を定義
        self.variable_values = self.define_variables()
        # ラベルと単位を定義
        self.variable_labels, self.unit_labels = self.define_labels()
        # h_dictの初期化
        self.h_dict = self.initialize_h_dict(x3, x4)
        # グリッド作成
        self.X, self.Y = self.create_grid()
        # 変数の値を代入
        self.h_values = self.assign_variable_values()
        # hの計算
        self.Z = self.calculate_h()

    # 変数を定義する関数
    def define_variables(self):
        return {
            "x1": np.linspace(0, 800, 200), # 溶接速度
            "x2": np.linspace(-6, 6, 200), # ヘッド位置
            "x3": np.linspace(0, 7, 200), # 材料の厚さ
            "x4": np.linspace(0, 2000, 200) # レーザーパワー
        }

    # ラベルと単位を定義する関数
    def define_labels(self):
        variable_labels = {
            "x1": "welding speed",
            "x2": "head position",
            "x3": "material thickness",
            "x4": "laser power"
        }

        unit_labels = {
            "x1": "mm/sec",
            "x2": "mm",
            "x3": "mm",
            "x4": "W"
        }

        return variable_labels, unit_labels

    # hの計算式
    def h(self, x1, x2, x3, x4):
        g_x1 =  np.exp(-self.lambda_1 * abs(x1))
        f_x2 = (1 / math.sqrt(2 * math.pi * self.sigma**2)) * np.exp(-((x2 - self.mu)**2) / (2 * self.sigma**2))
        s_x3 = self.C + (1 - self.C) * np.exp(-self.lambda_3 * x3)
        l_x4 = 1 + self.alpha * (x4 / 100)

        # return self.A1 * (f_x2 + g_x1 + s_x3 + l_x4)      # 線形結合
        return self.A1 * f_x2 * g_x1 * s_x3 * l_x4        # 乗法結合


        # h_dictの初期化
    def initialize_h_dict(self, x3, x4):
        return {
            "x1": self.variable_values["x1"], # 溶接速度
            "x2": self.variable_values["x2"], # ヘッド位置
            "x3": x3, # 材料の厚さ
            "x4": x4  # レーザーパワー
        }

    # グリッドを作成する関数
    def create_grid(self):
        X, Y = np.meshgrid(self.h_dict["x1"], self.h_dict["x2"])
        return X, Y

    # 変数の値を代入する関数
    def assign_variable_values(self):
        x1 = self.X
        x2 = self.Y
        x3 = self.h_dict["x3"]
        x4 = self.h_dict["x4"]
        h_values = [x1, x2, x3, x4]
        return h_values

    # hを計算する関数
    def calculate_h(self):
        return self.h(self.h_values[0], self.h_values[1], self.h_values[2], self.h_values[3])

    # フィギュアを取得する関数
    def get_figure(self):
        # Figure オブジェクトの作成
        figure = Figure(figsize=(6.4, 4.8), dpi=100)
        # サブプロットの追加
        ax = figure.add_subplot(111)

        # データの可視化
        cax = ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, vmax=self.amp)
        # カラーバーの追加
        figure.colorbar(cax, ax=ax, label='Amplitude')
        # グリッドの表示
        ax.grid()

        # x軸とy軸のラベル
        ax.set_xlabel(f"{self.variable_labels['x1']} [{self.unit_labels['x1']}]")
        ax.set_ylabel(f"{self.variable_labels['x2']} [{self.unit_labels['x2']}]")

        # タイトルの設定
        title_parts = []
        title_parts.append(f"{self.variable_labels['x3']}: {self.h_dict['x3']} [{self.unit_labels['x3']}]")
        title_parts.append(f"{self.variable_labels['x4']}: {self.h_dict['x4']} [{self.unit_labels['x4']}]")
        ax.set_title(", ".join(title_parts))

        return figure

class WeldingSimulationGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Welding Simulation")
        self.geometry("660x680")

        # タブコントロールの作成
        self.tab_control = ttk.Notebook(self)
        self.tab_control.pack(expand=1, fill='both')

        # Mainタブの作成
        self.tab_main = ttk.Frame(self.tab_control)
        self.tab_control.add(self.tab_main, text='main')

        # Settingタブの作成
        self.tab_setting = ttk.Frame(self.tab_control)
        self.tab_control.add(self.tab_setting, text='setting')

        # Mainタブのスピンボックス設定
        self.create_spinboxes()

        # Calculationボタンの設定
        self.calc_button = ttk.Button(self.tab_main, text="Calculation", command=self.calculate_and_plot)
        self.calc_button.pack(padx=10, pady=10)

        # グラフ表示用キャンバスの設定
        self.canvas_frame = ttk.Frame(self.tab_main)
        self.canvas_frame.pack(padx=10, pady=10)

        # Exitボタンの設定
        self.exit_button = ttk.Button(self.tab_main, text="EXIT", command=self.destroy)
        self.exit_button.pack(padx=10, pady=10)

        # Settingタブのエントリーボックス設定
        self.create_entry_boxes()

        # self.parametersの初期化
        self.parameters = {name: var.get() for name, var in self.setting_vars.items()}

        # 適応ボタンの設定
        self.apply_button = ttk.Button(self.tab_setting, text="Apply", command=self.apply_changes)
        self.apply_button.pack(padx=10, pady=10)

        # デフォルト値でグラフを表示
        self.calculate_and_plot()


    def create_spinboxes(self):
        # スピンボックスとラベルを作成
        self.var = [tk.DoubleVar(), tk.IntVar()]
        self.var[0].set(1.0)
        self.var[1].set(1000)

        self.labels = ["material thickness (mm)", "laser power (W)"]
        self.configs = [{"from_": 0, "to": 7.0, "increment": 0.1}, {"from_": 0, "to": 2000, "increment": 10}]

        # 共通のフレームを作成
        frame = ttk.Frame(self.tab_main)
        frame.pack(fill='x', padx=10, pady=5)

        for i, config in enumerate(self.configs):
            label = ttk.Label(frame, text=self.labels[i])
            # gridメソッドでラベルを配置
            label.grid(row=0, column=2*i, padx=10)
            spinbox = ttk.Spinbox(frame, from_=config["from_"], to=config["to"], increment=config["increment"], textvariable=self.var[i])
            # gridメソッドでスピンボックスを配置
            spinbox.grid(row=0, column=2*i+1, padx=10)

    def create_entry_boxes(self):
        # エントリーボックスとラベルを作成
        self.setting_vars = {"A1": tk.DoubleVar(value=1.0), 
                            "lambda_1": tk.DoubleVar(value=0.0045), 
                            "mu": tk.DoubleVar(value=0),
                            "sigma": tk.DoubleVar(value=2),
                            "lambda_3": tk.DoubleVar(value=1.0),
                            "C": tk.DoubleVar(value=0.5),
                            "alpha": tk.DoubleVar(value=1.0),
                            "amp": tk.DoubleVar(value=1.0)}

        for i, (label_text, var) in enumerate(self.setting_vars.items()):
            frame = ttk.Frame(self.tab_setting)
            frame.pack(fill='x', padx=10, pady=5)
            label = ttk.Label(frame, text=label_text, width=10)
            label.pack(side='left')
            entry_box = ttk.Entry(frame, textvariable=var)
            entry_box.pack(side='left')

    def apply_changes(self):
        self.parameters = {name: var.get() for name, var in self.setting_vars.items()}
        self.calculate_and_plot()

    def calculate_and_plot(self):
        x3 = self.var[0].get()
        x4 = self.var[1].get()

        A1 = self.parameters["A1"]
        lambda_1 = self.parameters["lambda_1"]
        mu = self.parameters["mu"]
        sigma = self.parameters["sigma"]
        lambda_3 = self.parameters["lambda_3"]
        C = self.parameters["C"]
        alpha = self.parameters["alpha"]
        amp = self.parameters["amp"]

        welding_simulation = WeldingSimulation(x3, x4, A1, lambda_1, mu, sigma, lambda_3, C, alpha, amp)
        figure = welding_simulation.get_figure()

        # 以前のグラフを削除
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()

        # グラフをキャンバスに表示
        self.canvas = FigureCanvasTkAgg(figure, self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # 新しいグラフを追加
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    app = WeldingSimulationGUI()
    app.mainloop()

