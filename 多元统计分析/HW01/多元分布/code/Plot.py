import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
 
# 支持中文字体
rcParams['font.sans-serif'] = ['SimHei']
rcParams['axes.unicode_minus'] = False

# 读取文件
file_path = 'Data/OXY.csv'
data = pd.read_csv(file_path)

# OXY vs time
plt.figure(figsize=(8, 6))
plt.scatter(data['time'], data['OXY'])
plt.title('OXY vs Time')
plt.xlabel('Time')
plt.ylabel('OXY')
plt.grid(True)
plt.savefig('figure/OXY_vs_Time.png')
plt.show()

# OXY vs age
plt.figure(figsize=(8, 6))
plt.scatter(data['age'], data['OXY'])
plt.title('OXY vs Age')
plt.xlabel('Age')
plt.ylabel('OXY')
plt.grid(True)
plt.savefig('figure/OXY_vs_Age.png')
plt.show()

# 选择指标
data_selected = data[['age', 'weight', 'time', 'spulse', 'rpulse', 'mpulse', 'OXY']]

# 绘制散点图矩阵
sns.pairplot(data_selected)
plt.suptitle("七项指标的散点图矩阵", y=1.02)
plt.savefig('figure/Scatter_Plot_Matrix.png')
plt.show()

# 获取序号为 1, 2, 21, 22 的数据
selected_samples = data.loc[data['序号'].isin([1, 2, 21, 22])]

# 获取相关特征
features = ['age', 'weight', 'time', 'spulse', 'rpulse', 'mpulse', 'OXY']
selected_data = selected_samples[features]

# 标准化
standardized_data = (selected_data - selected_data.mean()) / selected_data.std()

# 绘制轮廓图
plt.figure(figsize=(10, 6))
for i, sample in enumerate(standardized_data.values):
    plt.plot(features, sample, marker='o', label=f'样本 {selected_samples.iloc[i]["序号"]}')

plt.title('轮廓图')
plt.xlabel('7项指标')
plt.ylabel('标准化数据')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('figure/Contour_Plot.png')
plt.show()

# 雷达图绘制
labels = features
num_vars = len(labels)

angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

plt.figure(figsize=(8, 8))

for i, sample in enumerate(standardized_data.values):
    values = sample.tolist()
    values += values[:1]

    plt.polar(angles, values, marker='o', label=f'样本 {selected_samples.iloc[i]["序号"]}')

# 绘制雷达图
plt.title('雷达图')
plt.xticks(angles[:-1], labels, fontsize=12)
plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1))
plt.grid(True)
plt.tight_layout()
plt.savefig('figure/Radar_Chart.png')
plt.show()

# 调和曲线图函数
def Andrews(x):
    if isinstance(x, pd.DataFrame):
        x = x.values

    # 设置调和曲线自变量
    t = np.linspace(-np.pi, np.pi, 500)

    # 获得数据维度
    m, n = x.shape

    # 储存数据
    f = np.zeros((m, len(t)))

    # 代入调和曲线表达式
    for i in range(m):
        f[i, :] = x[i, 0] / np.sqrt(2)
        for j in range(1, n):
            if j % 2 == 0:
                f[i, :] += x[i, j] * np.cos(j // 2 * t)
            else:
                f[i, :] += x[i, j] * np.sin(j // 2 * t)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot([np.min(t), np.max(t)], [np.min(f), np.max(f)], 'k-', alpha=0)  # To set the plot limits
    plt.title('调和曲线图', fontsize=16)
    plt.xlabel('t', fontsize=14)
    plt.ylabel('f(t)', fontsize=14)

    # Plot each row in f as a separate line
    for i in range(m):
        plt.plot(t, f[i, :], label=f'样本 {i + 1}')

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('figure/Harmonic_Curve.png')
    plt.show()


Andrews(standardized_data)
