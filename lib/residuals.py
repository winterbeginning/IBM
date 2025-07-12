#!/usr/bin/env python3
import os
import matplotlib
matplotlib.rc("font", family='WenQuanYi Micro Hei')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import re
from collections import defaultdict
import matplotlib.cm as cm

# 配置参数
FILENAME = "postProcessing/residuals/0/residuals.dat"  # 数据文件名
REFRESH_INTERVAL = 2000     # 刷新间隔(毫秒)

class DataMonitor:
    def __init__(self, filename):
        self.filename = filename
        self.last_size = 0
        self.data = defaultdict(list)
        self.headers = []
        self.column_indices = {}
        self.color_map = {}
        self.lines = {}
        self.file_initialized = False
        
        # 初始化图形
        self.fig, self.ax = plt.subplots(figsize=(12, 6))
        plt.subplots_adjust(right=0.85)
        
        # 设置图形属性
        self.ax.set_xlabel('时间')
        self.ax.set_ylabel('残差')
        self.ax.set_title(f'实时监控: {os.path.basename(filename)}')
        self.ax.grid(True)
        self.ax.set_yscale('log')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(1e-6, 1)
        self.fig.canvas.manager.set_window_title('实时残差监控')
        
        # 状态文本
        self.status_text = self.ax.text(0.02, 0.98, '等待数据...', 
                                        transform=self.ax.transAxes,
                                        verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # 颜色映射
        self.colors = cm.tab10(np.linspace(0, 1, 10))
    
    def initialize_file(self):
        """初始化文件读取状态"""
        if os.path.exists(self.filename):
            try:
                self.last_size = 0
                self.read_all_data()
                self.file_initialized = True
                self.status_text.set_text('文件初始化完成')
                return True
            except Exception as e:
                self.status_text.set_text(f'初始化错误: {str(e)}')
                return False
        self.status_text.set_text('文件不存在')
        return False
    
    def read_data(self, all_lines=False):
        """读取文件数据"""
        if not os.path.exists(self.filename):
            return False
        
        try:
            with open(self.filename, 'r') as f:
                if all_lines:
                    self.last_size = 0
                    lines = f.readlines()
                else:
                    f.seek(self.last_size)
                    lines = f.readlines()
                self.last_size = f.tell()
        except (OSError, PermissionError) as e:
            self.status_text.set_text(f'读取错误: {str(e)}')
            return False
        
        if not lines:
            return False
        
        # 处理行数据
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                if 'Time' in line or 'time' in line.lower():
                    self.process_header(line)
                continue
            
            self.process_data_line(line)
        
        return True
    
    def read_all_data(self):
        """读取整个文件内容"""
        return self.read_data(all_lines=True)
    
    def read_new_data(self):
        """读取新增数据"""
        try:
            current_size = os.path.getsize(self.filename)
        except OSError as e:
            self.status_text.set_text(f'文件访问错误: {str(e)}')
            return False
        
        # 处理文件重置
        if current_size < self.last_size:
            self.last_size = 0
            for key in self.data:
                self.data[key] = []
            self.status_text.set_text('文件重置，重新开始读取...')
        
        # 检查文件变化
        if current_size == self.last_size:
            if not self.file_initialized:
                return self.initialize_file()
            return True
        
        return self.read_data()
    
    def process_header(self, header_line):
        """处理头部信息"""
        clean_line = re.sub(r'[^a-zA-Z0-9\s]', ' ', header_line[1:])
        headers = clean_line.split()
        
        if not headers:
            return
            
        # 确保时间列为第一列
        if headers[0].lower() not in ['time', 't']:
            headers.insert(0, 'Time')
        
        self.headers = headers
        self.column_indices = {name: idx for idx, name in enumerate(headers)}
        
        # 为列分配颜色
        for idx, col in enumerate(self.headers[1:]):
            self.color_map[col] = self.colors[idx % len(self.colors)]
        
        self.status_text.set_text(f'发现{len(self.headers)-1}个数据列')
    
    def process_data_line(self, line):
        """处理数据行"""
        parts = line.split()
        if not parts:
            return
        
        # 处理缺失头部的情况
        if not self.headers:
            try:
                float(parts[0])
                self.headers = [f'列_{i}' for i in range(len(parts))]
                self.headers[0] = '时间'
                self.column_indices = {name: idx for idx, name in enumerate(self.headers)}
                for idx, col in enumerate(self.headers[1:]):
                    self.color_map[col] = self.colors[idx % len(self.colors)]
                self.status_text.set_text(f'推断{len(self.headers)-1}个数据列')
            except ValueError:
                return
        
        # 填充缺失值
        if len(parts) < len(self.headers):
            parts.extend(['0'] * (len(self.headers) - len(parts)))
        
        # 处理时间数据
        try:
            time_val = float(parts[0])
            time_key = '时间' if '时间' in self.headers else 'Time'
            self.data[time_key].append(time_val)
        except ValueError:
            return
        
        # 处理其他列数据
        for col in self.headers[1:]:
            idx = self.column_indices[col]
            if idx < len(parts):
                val_str = parts[idx]
                try:
                    if val_str in ['N/A', 'NaN', 'nan', '-']:
                        self.data[col].append(np.nan)
                    else:
                        self.data[col].append(float(val_str))
                except ValueError:
                    self.data[col].append(np.nan)
    
    def update_plot(self, frame):
        """更新图表"""
        try:
            has_new_data = self.read_new_data()
            time_key = '时间' if '时间' in self.data else 'Time'
            
            if not self.data.get(time_key, []):
                return [self.status_text]
            
            # 更新曲线
            for col in self.headers[1:]:
                if col in self.data and self.data[col]:
                    min_len = min(len(self.data[time_key]), len(self.data[col]))
                    if min_len > 0:
                        x = self.data[time_key][:min_len]
                        y = self.data[col][:min_len]
                        
                        if col not in self.lines:
                            self.lines[col], = self.ax.plot(
                                x, y, '-', 
                                color=self.color_map.get(col, 'blue'),
                                alpha=0.8,
                                label=col
                            )
                        else:
                            self.lines[col].set_data(x, y)
            
            # 更新状态
            last_time = self.data[time_key][-1] if self.data[time_key] else 0
            status = f'时间: {last_time:.4e}\n数据点: {len(self.data[time_key])}'
            status += "\n文件已完整" if not has_new_data and self.file_initialized else "\n监控中..."
            self.status_text.set_text(status)
            
            # 调整坐标轴
            if self.data[time_key]:
                x_min, x_max = min(self.data[time_key]), max(self.data[time_key])
                x_range = max(x_max - x_min, 1e-10, x_max * 0.1)
                self.ax.set_xlim(x_min - 0.05 * x_range, x_max + 0.05 * x_range)
                
                # 计算Y轴范围
                y_vals = []
                for col in self.headers[1:]:
                    if col in self.data:
                        y_vals.extend(y for y in self.data[col] if not np.isnan(y) and y > 0)
                
                if y_vals:
                    y_min, y_max = min(y_vals), max(y_vals)
                    y_min_display = max(1e-10, y_min * 0.1)
                    y_max_display = y_max * 10
                    if y_min_display < y_max_display:
                        self.ax.set_ylim(y_min_display, y_max_display)
            
            # 更新图例
            if self.lines:
                handles = list(self.lines.values())[:10]
                self.ax.legend(handles, [line.get_label() for line in handles], 
                              loc='upper right', fontsize=9)
            
            return [self.status_text] + list(self.lines.values())
        
        except Exception as e:
            self.status_text.set_text(f'错误: {str(e)}')
            return [self.status_text]
    
    def start_monitoring(self):
        """启动监控"""
        self.initialize_file()
        # 创建动画对象但不调用start()
        self.ani = FuncAnimation(
            self.fig, 
            self.update_plot, 
            interval=REFRESH_INTERVAL,
            blit=True,
            cache_frame_data=False
        )
        plt.show()  # 直接调用plt.show()启动动画

def main():
    print("启动实时残差监控...")
    print(f"监控文件: {FILENAME}")
    print(f"刷新间隔: {REFRESH_INTERVAL}毫秒")
    print("按Ctrl+C退出")
    
    try:
        DataMonitor(FILENAME).start_monitoring()
    except KeyboardInterrupt:
        print("\n监控已停止")
    except Exception as e:
        print(f"发生错误: {str(e)}")

if __name__ == "__main__":
    main()