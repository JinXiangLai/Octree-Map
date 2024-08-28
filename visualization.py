import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# point_cloud_file = './build/point_cloud_full_octree.csv'
point_cloud_file = './build/point_cloud_dynamic_create_octree.csv'
# point_cloud_file = './build/point_cloud_2_make_octree.csv'

# points是一个Nx3的数组，其中N是点的数量，3代表x, y, z坐标
points = []
with open(point_cloud_file, 'r') as file:
    for line in file:
        x, y, z = map(float, line.split(' '))
        points.append([x, y, z])
points = np.array(points)
 
# 显示点云大小
point_size = 0.1
# 计算中心点
center = points.mean(axis=0)
print('center: ', center)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:, 0], points[:, 1], points[:, 2], s = point_size)

# 设置视图，‌使其以中心点为中心
lim = 5
ax.set_xlim(center[0] - lim, center[0] + lim)
ax.set_ylim(center[1] - lim, center[1] + lim)
ax.set_zlim(center[2] - lim, center[2] + lim)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()