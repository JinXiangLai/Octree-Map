import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

VISUALIZE_SCAN2MAP = 1

def get_point_cloud(point_cloud_file:str) -> np.ndarray:
    # points是一个Nx3的数组，其中N是点的数量，3代表x, y, z坐标
    points = []
    with open(point_cloud_file, 'r') as file:
        for line in file:
            x, y, z = map(float, line.split(' '))
            points.append([x, y, z])
    return np.array(points)

def show(ax):
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

if not VISUALIZE_SCAN2MAP:
    # point_cloud_file = './build/point_cloud_full_octree.csv'
    # point_cloud_file = './build/point_cloud_2_make_octree.csv'
    point_cloud_file = './build/point_cloud_dynamic_create_octree.csv'
    
    points = get_point_cloud(point_cloud_file)

    # 显示点云大小
    point_size = 0.1
    # 计算中心点
    center = points.mean(axis=0)
    print('center: ', center)
    # 配置点颜色
    colors = np.random.rand(points.shape[0], 3)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s = point_size, c = colors)

    # 设置视图，‌使其以中心点为中心
    lim = 500
    ax.set_xlim(center[0] - lim, center[0] + lim)
    ax.set_ylim(center[1] - lim, center[1] + lim)
    ax.set_zlim(center[2] - lim, center[2] + lim)

    show(ax)

else:
    scan_cloud_file = './build/scan_cloud.csv'
    local_map_file = './build/local_map.csv'
    scan_point = get_point_cloud(scan_cloud_file)
    map_point = get_point_cloud(local_map_file)
    scan_size = 1
    scan_color = 'r'
    map_size = 0.1
    map_color = 'b'
    center = map_point.mean(axis=0)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(map_point[:, 0], map_point[:, 1], map_point[:, 2], s = map_size, c = map_color)
    ax.scatter(scan_point[:, 0], scan_point[:, 1], scan_point[:, 2], s = scan_size, c = scan_color)
    show(ax)
