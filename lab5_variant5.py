import tkinter as tk
from tkinter import filedialog, messagebox

class ClippingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Лабораторная работа 5 - Вариант 5")
        self.root.geometry("1000x700")
        
        # Данные
        self.segments = []
        self.polygons = []
        self.clipping_window = None  # Прямоугольное окно отсечения
        self.clipped_segments = []
        self.clipped_polygons = []
        self.clipping_mode = tk.StringVar(value="segments")  # "segments", "polygons"
        self.show_without_clipping = tk.BooleanVar(value=False)  # Показать без отсечения
        
        # Параметры масштабирования для координатной системы
        self.scale_params = {
            'scale': 1.0,
            'center_x': 0.0,
            'center_y': 0.0,
            'origin_screen_x': 0.0,
            'origin_screen_y': 0.0
        }
        
        # Создание интерфейса
        self.create_ui()
        
    def create_ui(self):
        # Панель управления
        control_frame = tk.Frame(self.root)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        
        tk.Button(control_frame, text="Загрузить файл", command=self.load_file).pack(side=tk.LEFT, padx=5)
        tk.Button(control_frame, text="Очистить", command=self.clear).pack(side=tk.LEFT, padx=5)
        
        # Выбор метода отсечения
        method_frame = tk.LabelFrame(control_frame, text="Метод отсечения", padx=5, pady=5)
        method_frame.pack(side=tk.LEFT, padx=10)
        
        tk.Radiobutton(method_frame, text="Отрезки (Лианга-Барски)", 
                      variable=self.clipping_mode, value="segments",
                      command=self.on_mode_change).pack(anchor=tk.W)
        tk.Radiobutton(method_frame, text="Многоугольники", 
                      variable=self.clipping_mode, value="polygons",
                      command=self.on_mode_change).pack(anchor=tk.W)
        
        # Чекбокс для отображения без отсечения
        tk.Checkbutton(control_frame, text="Показать без отсечения", 
                      variable=self.show_without_clipping,
                      command=self.on_mode_change).pack(side=tk.LEFT, padx=10)
        
        # Canvas для отрисовки
        self.canvas = tk.Canvas(self.root, bg="white", width=950, height=600)
        self.canvas.pack(padx=5, pady=5)
        
        # Информационная панель
        self.info_label = tk.Label(self.root, text="Загрузите файл с данными", anchor=tk.W)
        self.info_label.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=5)
    
    def on_mode_change(self):
        """Вызывается при изменении режима отсечения"""
        if self.clipping_window:
            self.draw()
        
    def load_file(self):
        filename = filedialog.askopenfilename(
            title="Выберите файл с данными",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            try:
                self.parse_file(filename)
                self.draw()
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при загрузке файла: {str(e)}")
    
    def parse_file(self, filename):
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        if len(lines) < 2:
            raise ValueError("Файл должен содержать минимум 2 строки")
        
        n = int(lines[0])
        
        if len(lines) < n + 2:
            raise ValueError(f"Недостаточно строк в файле. Ожидается {n + 2}, получено {len(lines)}")
        
        # Читаем отрезки/многоугольники
        self.segments = []
        self.polygons = []
        
        for i in range(1, n + 1):
            coords = list(map(float, lines[i].split()))
            if len(coords) == 4:
                # Отрезок (4 координаты: x1 y1 x2 y2)
                self.segments.append({
                    'x1': coords[0], 'y1': coords[1],
                    'x2': coords[2], 'y2': coords[3]
                })
            elif len(coords) >= 6 and len(coords) % 2 == 0:
                # Многоугольник (четное количество координат, минимум 6 для треугольника)
                polygon = []
                for j in range(0, len(coords), 2):
                    polygon.append((coords[j], coords[j+1]))
                self.polygons.append(polygon)
            else:
                raise ValueError(f"Неверный формат координат в строке {i+1}: ожидается 4 координаты для отрезка или четное количество >= 6 для многоугольника")
        
        # Читаем отсекающее окно (прямоугольное)
        window_coords = list(map(float, lines[n + 1].split()))
        
        if len(window_coords) != 4:
            raise ValueError("Отсекающее окно должно содержать 4 координаты: xmin ymin xmax ymax")
        
        self.clipping_window = {
            'xmin': window_coords[0],
            'ymin': window_coords[1],
            'xmax': window_coords[2],
            'ymax': window_coords[3]
        }
        
        # Проверяем корректность окна
        if self.clipping_window['xmin'] >= self.clipping_window['xmax'] or \
           self.clipping_window['ymin'] >= self.clipping_window['ymax']:
            raise ValueError("Неверные координаты отсекающего окна: xmin < xmax и ymin < ymax")
        
        # Вычисляем параметры масштабирования
        self.calculate_scale_params()
        
        # Выполняем отсечение
        self.perform_clipping()
        
        # Формируем информационное сообщение
        cw = self.clipping_window
        window_info = f"Прямоугольное окно: ({cw['xmin']}, {cw['ymin']}) - ({cw['xmax']}, {cw['ymax']})"
        
        self.info_label.config(text=f"Загружено: {len(self.segments)} отрезков, {len(self.polygons)} многоугольников. {window_info}")
    
    def liang_barsky_clip(self, x1, y1, x2, y2, xmin, ymin, xmax, ymax):
        """Алгоритм Лианга-Барски для отсечения отрезка"""
        dx = x2 - x1
        dy = y2 - y1
        
        p = [-dx, dx, -dy, dy]
        q = [x1 - xmin, xmax - x1, y1 - ymin, ymax - y1]
        
        u1 = 0.0
        u2 = 1.0
        
        for i in range(4):
            if p[i] == 0:
                # Линия параллельна границе
                if q[i] < 0:
                    return None  # Отрезок полностью вне окна
            else:
                r = q[i] / p[i]
                if p[i] < 0:
                    u1 = max(u1, r)
                else:
                    u2 = min(u2, r)
        
        if u1 > u2:
            return None  # Отрезок полностью вне окна
        
        # Вычисляем конечные точки отсеченного отрезка
        if u1 > 0 or u2 < 1:
            x1_clipped = x1 + u1 * dx
            y1_clipped = y1 + u1 * dy
            x2_clipped = x1 + u2 * dx
            y2_clipped = y1 + u2 * dy
            return (x1_clipped, y1_clipped, x2_clipped, y2_clipped)
        
        # Отрезок полностью внутри
        return (x1, y1, x2, y2)
    
    def polygon_area(self, polygon):
        """Возвращает удвоенную площадь многоугольника (знак определяет ориентацию)."""
        area = 0.0
        n = len(polygon)
        for i in range(n):
            x1, y1 = polygon[i]
            x2, y2 = polygon[(i + 1) % n]
            area += x1 * y2 - x2 * y1
        return area / 2.0

    def ensure_ccw(self, polygon):
        """Гарантирует обход вершин против часовой стрелки."""
        if self.polygon_area(polygon) < 0:
            return list(reversed(polygon))
        return polygon

    def is_convex(self, polygon):
        """Проверяет, является ли многоугольник выпуклым"""
        if len(polygon) < 3:
            return False
        
        n = len(polygon)
        sign = 0
        
        for i in range(n):
            p1 = polygon[i]
            p2 = polygon[(i + 1) % n]
            p3 = polygon[(i + 2) % n]
            
            # Вычисляем векторное произведение
            cross = (p2[0] - p1[0]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[0] - p2[0])
            
            if cross != 0:
                if sign == 0:
                    sign = 1 if cross > 0 else -1
                elif (sign > 0 and cross < 0) or (sign < 0 and cross > 0):
                    return False
        
        return True
    
    def point_in_polygon(self, point, polygon):
        """Проверяет, находится ли точка внутри выпуклого многоугольника (полученного в CCW)."""
        if len(polygon) < 3:
            return False
        
        poly_ccw = self.ensure_ccw(polygon)
        x, y = point
        n = len(poly_ccw)
        for i in range(n):
            x1, y1 = poly_ccw[i]
            x2, y2 = poly_ccw[(i + 1) % n]
            # cross >= 0 для точки слева от ребра (внутри для CCW)
            cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
            if cross < -1e-10:
                return False
        return True
    
    def clip_polygon(self, polygon, xmin, ymin, xmax, ymax):
        """Алгоритм отсечения выпуклого многоугольника (Sutherland-Hodgman)"""
        # Проверяем выпуклость
        if not self.is_convex(polygon):
            # Для невыпуклых многоугольников возвращаем пустой результат
            # или можно использовать другой алгоритм
            return []
        def inside(point, edge):
            x, y = point
            if edge == 'left':
                return x >= xmin
            elif edge == 'right':
                return x <= xmax
            elif edge == 'bottom':
                return y >= ymin
            elif edge == 'top':
                return y <= ymax
        
        def intersect(p1, p2, edge):
            x1, y1 = p1
            x2, y2 = p2
            
            if edge == 'left':
                if x1 == x2:
                    return None
                t = (xmin - x1) / (x2 - x1)
                y = y1 + t * (y2 - y1)
                if ymin <= y <= ymax:
                    return (xmin, y)
            elif edge == 'right':
                if x1 == x2:
                    return None
                t = (xmax - x1) / (x2 - x1)
                y = y1 + t * (y2 - y1)
                if ymin <= y <= ymax:
                    return (xmax, y)
            elif edge == 'bottom':
                if y1 == y2:
                    return None
                t = (ymin - y1) / (y2 - y1)
                x = x1 + t * (x2 - x1)
                if xmin <= x <= xmax:
                    return (x, ymin)
            elif edge == 'top':
                if y1 == y2:
                    return None
                t = (ymax - y1) / (y2 - y1)
                x = x1 + t * (x2 - x1)
                if xmin <= x <= xmax:
                    return (x, ymax)
            return None
        
        # Применяем отсечение по каждой границе
        edges = ['left', 'right', 'bottom', 'top']
        result = polygon.copy()
        
        for edge in edges:
            if not result or len(result) < 3:
                break
            new_polygon = []
            n = len(result)
            
            # Берем последнюю точку как предыдущую
            prev_point = result[-1]
            prev_inside = inside(prev_point, edge)
            
            for i in range(n):
                current = result[i]
                current_inside = inside(current, edge)
                
                if prev_inside and current_inside:
                    # Обе точки внутри - добавляем текущую
                    new_polygon.append(current)
                elif prev_inside and not current_inside:
                    # Предыдущая внутри, текущая снаружи - добавляем пересечение
                    intersection = intersect(prev_point, current, edge)
                    if intersection:
                        new_polygon.append(intersection)
                elif not prev_inside and current_inside:
                    # Предыдущая снаружи, текущая внутри - добавляем пересечение и текущую
                    intersection = intersect(prev_point, current, edge)
                    if intersection:
                        new_polygon.append(intersection)
                    new_polygon.append(current)
                # Если обе снаружи - ничего не добавляем
                
                prev_point = current
                prev_inside = current_inside
            
            result = new_polygon
        
        return result if len(result) >= 3 else []
    
    def perform_clipping(self):
        """Выполняет отсечение всех отрезков и многоугольников"""
        if not self.clipping_window:
            return
        
        cw = self.clipping_window
        
        # Отсечение отрезков алгоритмом Лианга-Барски
        self.clipped_segments = []
        for seg in self.segments:
            clipped = self.liang_barsky_clip(
                seg['x1'], seg['y1'], seg['x2'], seg['y2'],
                cw['xmin'], cw['ymin'], cw['xmax'], cw['ymax']
            )
            if clipped:
                self.clipped_segments.append(clipped)
        
        # Отсечение многоугольников алгоритмом Сазерленда-Ходжмана
        self.clipped_polygons = []
        non_convex_count = 0
        for poly in self.polygons:
            if not self.is_convex(poly):
                non_convex_count += 1
                continue  # Пропускаем невыпуклые многоугольники
            
            clipped = self.clip_polygon(
                poly, cw['xmin'], cw['ymin'], cw['xmax'], cw['ymax']
            )
            
            if clipped:
                self.clipped_polygons.append(clipped)
        
        if non_convex_count > 0:
            self.info_label.config(
                text=f"{self.info_label.cget('text')} | ВНИМАНИЕ: {non_convex_count} невыпуклых многоугольников пропущено"
            )
    
    def calculate_scale_params(self):
        """Вычисляет параметры масштабирования и сохраняет их"""
        # Находим границы всех объектов
        all_x = []
        all_y = []
        
        if self.clipping_window:
            cw = self.clipping_window
            all_x.extend([cw['xmin'], cw['xmax']])
            all_y.extend([cw['ymin'], cw['ymax']])
        
        for seg in self.segments:
            all_x.extend([seg['x1'], seg['x2']])
            all_y.extend([seg['y1'], seg['y2']])
        
        for poly in self.polygons:
            for px, py in poly:
                all_x.append(px)
                all_y.append(py)
        
        min_x, max_x = min(all_x), max(all_x)
        min_y, max_y = min(all_y), max(all_y)
        
        # Добавляем отступы
        padding = 50
        width = self.canvas.winfo_width() or 950
        height = self.canvas.winfo_height() or 600
        
        # Обработка случая, когда все точки совпадают
        range_x = max_x - min_x if max_x != min_x else max(abs(min_x), 1) * 2
        range_y = max_y - min_y if max_y != min_y else max(abs(min_y), 1) * 2
        
        # Если диапазон слишком мал, увеличиваем его
        if range_x < 0.1:
            range_x = 1
        if range_y < 0.1:
            range_y = 1
        
        scale_x = (width - 2 * padding) / range_x
        scale_y = (height - 2 * padding) / range_y
        scale = min(scale_x, scale_y)
        
        # Центрируем
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        
        origin_screen_x = width / 2 - center_x * scale
        origin_screen_y = height / 2 + center_y * scale
        
        # Сохраняем параметры
        self.scale_params = {
            'scale': scale,
            'center_x': center_x,
            'center_y': center_y,
            'origin_screen_x': origin_screen_x,
            'origin_screen_y': origin_screen_y,
            'width': width,
            'height': height
        }
    
    def world_to_screen(self, x, y):
        """Преобразует мировые координаты в экранные"""
        if not self.scale_params:
            return x, y
        
        sp = self.scale_params
        screen_x = sp['origin_screen_x'] + x * sp['scale']
        screen_y = sp['origin_screen_y'] - y * sp['scale']  # Инвертируем Y
        
        return screen_x, screen_y
    
    def draw(self):
        """Отрисовывает все объекты"""
        self.canvas.delete("all")
        
        if not self.clipping_window:
            # Если нет окна, все равно рисуем объекты
            if not self.segments and not self.polygons:
                return
        
        mode = self.clipping_mode.get()
        show_clipping = not self.show_without_clipping.get()
        
        # Рисуем координатные оси
        self.draw_axes()
        
        # Рисуем отсекающее окно (синий)
        if show_clipping and self.clipping_window:
            cw = self.clipping_window
            x1, y1 = self.world_to_screen(cw['xmin'], cw['ymin'])
            x2, y2 = self.world_to_screen(cw['xmax'], cw['ymax'])
            try:
                self.canvas.create_rectangle(x1, y1, x2, y2, outline="blue", width=2, 
                                            fill="lightblue", stipple="gray25")
            except:
                self.canvas.create_rectangle(x1, y1, x2, y2, outline="blue", width=2, fill="lightblue")
        
        # Рисуем исходные многоугольники (красный пунктир)
        if mode == "polygons":
            for poly in self.polygons:
                screen_points = [self.world_to_screen(px, py) for px, py in poly]
                if len(screen_points) >= 3:
                    # Рисуем исходный многоугольник
                    self.canvas.create_polygon(screen_points, outline="red", fill="", width=2, dash=(5, 5))
                    # Показываем вершины исходного многоугольника
                    for px, py in poly:
                        sx, sy = self.world_to_screen(px, py)
                        # Проверяем, внутри ли окна (если окно задано)
                        inside = False
                        if show_clipping and self.clipping_window:
                            cw = self.clipping_window
                            inside = (cw['xmin'] <= px <= cw['xmax'] and cw['ymin'] <= py <= cw['ymax'])
                        
                        if inside:
                            self.canvas.create_oval(sx - 4, sy - 4, sx + 4, sy + 4, 
                                                  fill="orange", outline="red", width=1)
                        else:
                            self.canvas.create_oval(sx - 3, sy - 3, sx + 3, sy + 3, 
                                                  fill="red", outline="darkred", width=1)
        
        # Рисуем отсеченные многоугольники (зеленый) - алгоритм отсечения многоугольника
        if mode == "polygons" and show_clipping:
            for i, poly in enumerate(self.clipped_polygons):
                screen_points = [self.world_to_screen(px, py) for px, py in poly]
                if len(screen_points) >= 3:
                    # Закрашиваем отсеченную часть ярким цветом для наглядности
                    # Это видимая часть после отсечения
                    self.canvas.create_polygon(screen_points, outline="darkgreen", fill="lightgreen", width=3)
                    # Показываем вершины отсеченного многоугольника
                    for px, py in poly:
                        sx, sy = self.world_to_screen(px, py)
                        self.canvas.create_oval(sx - 5, sy - 5, sx + 5, sy + 5, 
                                              fill="green", outline="darkgreen", width=2)
        
        # Рисуем исходные отрезки (красный пунктир)
        if mode == "segments":
            for seg in self.segments:
                x1, y1 = self.world_to_screen(seg['x1'], seg['y1'])
                x2, y2 = self.world_to_screen(seg['x2'], seg['y2'])
                self.canvas.create_line(x1, y1, x2, y2, fill="red", width=1, dash=(5, 5))
        
        # Рисуем отсеченные отрезки (зеленый жирный) - алгоритм Лианга-Барски
        if mode == "segments" and show_clipping:
            for clipped in self.clipped_segments:
                x1, y1 = self.world_to_screen(clipped[0], clipped[1])
                x2, y2 = self.world_to_screen(clipped[2], clipped[3])
                self.canvas.create_line(x1, y1, x2, y2, fill="green", width=3)
        
        # Рисуем легенду
        self.draw_legend()
    
    def draw_axes(self):
        """Рисует координатные оси с сеткой и подписями"""
        if not self.scale_params:
            return
        
        # Пересчитываем параметры масштабирования
        self.calculate_scale_params()
        
        if not self.scale_params:
            return
        
        sp = self.scale_params
        width = sp['width']
        height = sp['height']
        origin_screen_x = sp['origin_screen_x']
        origin_screen_y = sp['origin_screen_y']
        scale = sp['scale']
        
        # Находим границы для сетки
        all_x = []
        all_y = []
        
        if self.clipping_window:
            cw = self.clipping_window
            all_x.extend([cw['xmin'], cw['xmax']])
            all_y.extend([cw['ymin'], cw['ymax']])
        
        for seg in self.segments:
            all_x.extend([seg['x1'], seg['x2']])
            all_y.extend([seg['y1'], seg['y2']])
        
        for poly in self.polygons:
            for px, py in poly:
                all_x.append(px)
                all_y.append(py)
        
        if not all_x or not all_y:
            return
        
        min_x, max_x = min(all_x), max(all_x)
        min_y, max_y = min(all_y), max(all_y)
        
        range_x = max_x - min_x if max_x != min_x else 1
        range_y = max_y - min_y if max_y != min_y else 1
        
        if range_x < 0.1:
            range_x = 1
        if range_y < 0.1:
            range_y = 1
        
        # Вычисляем шаг сетки
        step = self.calculate_grid_step(range_x, range_y)
        
        # Рисуем сетку
        self.draw_grid(origin_screen_x, origin_screen_y, scale, step, width, height, min_x, max_x, min_y, max_y)
        
        # Рисуем оси
        self.canvas.create_line(0, origin_screen_y, width, origin_screen_y, fill="black", width=2)
        self.canvas.create_line(origin_screen_x, 0, origin_screen_x, height, fill="black", width=2)
        
        # Стрелки на осях
        arrow_size = 8
        # Стрелка оси X
        self.canvas.create_line(width - arrow_size, origin_screen_y - arrow_size, 
                               width, origin_screen_y, 
                               width - arrow_size, origin_screen_y + arrow_size, 
                               fill="black", width=2)
        # Стрелка оси Y
        self.canvas.create_line(origin_screen_x - arrow_size, arrow_size,
                               origin_screen_x, 0,
                               origin_screen_x + arrow_size, arrow_size,
                               fill="black", width=2)
        
        # Подписи осей
        self.canvas.create_text(width - 15, origin_screen_y - 15, text="X", 
                               fill="black", font=("Arial", 12, "bold"))
        self.canvas.create_text(origin_screen_x + 15, 15, text="Y", 
                               fill="black", font=("Arial", 12, "bold"))
        
        # Подпись начала координат
        if abs(origin_screen_x) < width and abs(origin_screen_y) < height:
            self.canvas.create_text(origin_screen_x - 10, origin_screen_y + 15, text="O", 
                                   fill="black", font=("Arial", 10, "bold"))
    
    def calculate_grid_step(self, range_x, range_y):
        """Вычисляет оптимальный шаг сетки"""
        max_range = max(range_x, range_y)
        
        # Выбираем шаг из стандартных значений
        if max_range <= 1:
            return 0.1
        elif max_range <= 5:
            return 0.5
        elif max_range <= 10:
            return 1
        elif max_range <= 50:
            return 5
        elif max_range <= 100:
            return 10
        elif max_range <= 500:
            return 50
        else:
            return 100
    
    def draw_grid(self, origin_x, origin_y, scale, step, width, height, min_x, max_x, min_y, max_y):
        """Рисует координатную сетку с подписями"""
        # Вертикальные линии сетки
        x = step * int(min_x / step)
        while x <= max_x:
            screen_x = origin_x + x * scale
            if 0 < screen_x < width:
                self.canvas.create_line(screen_x, 0, screen_x, height, 
                                       fill="lightgray", width=1, dash=(2, 2))
                # Подпись значения
                if abs(screen_x - origin_x) > 5:  # Не подписываем на оси
                    self.canvas.create_text(screen_x, origin_y + 15, text=f"{x:.1f}", 
                                           fill="gray", font=("Arial", 8))
            x += step
        
        # Горизонтальные линии сетки
        y = step * int(min_y / step)
        while y <= max_y:
            screen_y = origin_y - y * scale
            if 0 < screen_y < height:
                self.canvas.create_line(0, screen_y, width, screen_y, 
                                       fill="lightgray", width=1, dash=(2, 2))
                # Подпись значения
                if abs(screen_y - origin_y) > 5:  # Не подписываем на оси
                    self.canvas.create_text(origin_x - 10, screen_y, text=f"{y:.1f}", 
                                           fill="gray", font=("Arial", 8), anchor="e")
            y += step
    
    def draw_legend(self):
        """Рисует легенду"""
        x = 20
        y = 20
        
        mode = self.clipping_mode.get()
        
        # Отсекающее окно
        try:
            self.canvas.create_rectangle(x, y, x + 20, y + 15, outline="blue", fill="lightblue", stipple="gray25")
        except:
            self.canvas.create_rectangle(x, y, x + 20, y + 15, outline="blue", fill="lightblue")
        self.canvas.create_text(x + 30, y + 7, text="Отсекающее окно", anchor="w", font=("Arial", 9))
        
        # Исходные объекты
        y += 25
        if mode == "segments":
            self.canvas.create_line(x, y + 7, x + 20, y + 7, fill="red", width=1, dash=(5, 5))
            self.canvas.create_text(x + 30, y + 7, text="Исходные отрезки", anchor="w", font=("Arial", 9))
        elif mode == "polygons":
            self.canvas.create_polygon([x, y + 2, x + 20, y + 2, x + 10, y + 12], 
                                      outline="red", fill="", width=2, dash=(5, 5))
            self.canvas.create_text(x + 30, y + 7, text="Исходные многоугольники (оранжевые точки - внутри окна)", anchor="w", font=("Arial", 8))
        
        # Отсеченные объекты
        y += 25
        if mode == "segments":
            self.canvas.create_line(x, y + 7, x + 20, y + 7, fill="green", width=3)
            self.canvas.create_text(x + 30, y + 7, text="Отсеченные отрезки (Лианга-Барски)", anchor="w", font=("Arial", 9))
        elif mode == "polygons":
            self.canvas.create_polygon([x, y + 2, x + 20, y + 2, x + 10, y + 12], 
                                      outline="darkgreen", fill="lightgreen", width=3)
            self.canvas.create_text(x + 30, y + 7, text="Отсеченные многоугольники (зеленая заливка)", anchor="w", font=("Arial", 8))
    
    def clear(self):
        """Очищает canvas"""
        self.canvas.delete("all")
        self.segments = []
        self.polygons = []
        self.clipping_window = None
        self.clipped_segments = []
        self.clipped_polygons = []
        self.info_label.config(text="Очищено")

if __name__ == "__main__":
    root = tk.Tk()
    app = ClippingApp(root)
    root.mainloop()

