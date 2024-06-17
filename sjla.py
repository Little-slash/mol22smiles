import numpy as np

# 给定数据
V = 1.0  # 体积 (m^3)
T = 300.0  # 温度 (K)
k_B = 1.380649e-23  # 玻尔兹曼常数 (J/K)
tau = 1  # 特征时间 (s)

# 转换为 Pa
rmsd_xy = 67.6086 * 1e5
rmsd_xz = 67.2 * 1e5
rmsd_yz = 77.2579 * 1e5

# 计算初始自相关函数值
C_xy_0 = rmsd_xy**2
C_xz_0 = rmsd_xz**2
C_yz_0 = rmsd_yz**2

# 近似积分
eta_xy = (V / (k_B * T)) * C_xy_0 * tau
eta_xz = (V / (k_B * T)) * C_xz_0 * tau
eta_yz = (V / (k_B * T)) * C_yz_0 * tau

print("粗略估算粘度 η_xy =", eta_xy, "Pa·s")
print("粗略估算粘度 η_xz =", eta_xz, "Pa·s")
print("粗略估算粘度 η_yz =", eta_yz, "Pa·s")
