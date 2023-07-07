
import tkinter as tk
from tkinter import ttk
import sys
import 手


start_prob = sys.argv[1] if len(sys.argv) >= 2 else None
pixel_size = 1 # negative = division


root = tk.Tk()
root.geometry("1000x1000")


shortcuts = {}

def key(k):
    def d(f):
        global shortcuts
        shortcuts[k] = f
        return f
    return d


def on_press(event) -> None:
    if event.keysym in shortcuts:
        shortcuts[event.keysym]()
    else:
        print(f'pressed {repr(event.keysym)}')

root.bind('<Key>', on_press)


prob_ids = 手.all_prob_ids()


prob_var = tk.StringVar()
prob_cb = ttk.Combobox(root, textvariable=prob_var)
prob_cb['values'] = prob_ids
if start_prob is not None:
    prob_cb.set(start_prob)
else:
    prob_cb.current(0)
prob_cb['state'] = 'readonly'
prob_cb.pack()


# sol_ids = 手.get_all_sol_ids()
cur_sol_ids = []

sol_var = tk.StringVar()
sol_cb = ttk.Combobox(root, textvariable=sol_var)
sol_cb['values'] = []
sol_cb['state'] = 'readonly'
sol_cb.pack()


top_label = tk.Label(root)
top_label.pack()


canvas_w = 1000
canvas_h = 1000
canvas = tk.Canvas(root, width=canvas_w, height=canvas_h)
canvas.pack(fill=tk.BOTH)


coord_label = tk.Label(root)
coord_label.place(x=0, y=0)


score_label = tk.Label(root)
score_label.place(x=0, y=25)


def set_score(score) -> None:
    if score is None:
        score_label.config(text=f'score: -')
    else:
        score_label.config(text=f'score: {score}')


def tk_color_from_rgba(rgba) -> str:
    (r, g, b, _) = rgba
    return f'#{r:02x}{g:02x}{b:02x}'


def project_coord(x: float, y: float) -> tuple[int, int]:
    if pixel_size > 0:
        return (round(x*pixel_size), canvas_h - round(y*pixel_size))
    if pixel_size < 0:
        return (round(x/-pixel_size), canvas_h - round(y/-pixel_size))
    assert(False)


def on_mouse_move(event) -> None:
    (x, y) = project_coord(event.x, event.y)
    coord_label.config(text=f'x: {x}, y: {y}')

canvas.bind("<Motion>", on_mouse_move)


def fit_pixel_size(w: float, h: float) -> int:
    size = 1
    while size*w <= canvas_w and size*h <= canvas_h:
        size += 1
    while size > 1 and (size*w > canvas_w or size*h > canvas_h):
        size -= 1
    if size == 1 and (size*w > canvas_w or size*h > canvas_h):
        while w/size > canvas_w or h/size > canvas_h:
            size += 1
        return -size
    return size


def draw_prob(prob) -> None:
    global pixel_size
    canvas.delete('prob')

    rw = prob['room_width']
    rh = prob['room_height']
    pixel_size = fit_pixel_size(rw, rh)
    print(pixel_size)

    # room
    (x0, y0) = project_coord(0, 0)
    (x1, y1) = project_coord(rw, rh)
    canvas.create_rectangle(
        x0, y0, x1, y1,
        fill='white', width=0, tags='prob'
    )

    # stage
    (sx, sy) = prob['stage_bottom_left']
    (x0, y0) = project_coord(sx, sy)
    sw = prob['stage_width']
    sh = prob['stage_height']
    (x1, y1) = project_coord(sx+sw, sy+sh)
    canvas.create_rectangle(
        x0, y0, x1, y1,
        fill='white', outline='red', width=1, tags='prob'
    )

    # attendees
    for a in prob['attendees']:
        x = a['x']
        y = a['y']
        (x0, y0) = project_coord(x, y)
        (x1, y1) = project_coord(x, y)
        canvas.create_rectangle(
            x0, y0, x1, y1,
            fill='black', width=0, tags='prob'
        )


def update_top_label(prob) -> None:
    p = prob
    text = (
        f'musicians: {len(p["musicians"])}, '
        f'instuments: {len(set(p["musicians"]))}, '
        f'attendees: {len(p["attendees"])}, '
        f'room_size: {p["room_width"]}x{p["room_height"]}, '
        f'stage_size: {p["stage_width"]}x{p["stage_height"]}'
    )
    top_label.config(text=text)


def switch_prob(prob_id: int) -> None:
    prob = 手.get_prob(prob_id)
    draw_prob(prob)
    update_top_label(prob)
    coord_label.focus_set()


def switch_sol(sol_id) -> None:
    coord_label.focus_set()


def prob_changed(event) -> None:
    switch_prob(int(prob_cb.get()))

prob_cb.bind('<<ComboboxSelected>>', prob_changed)


def sol_changed(event) -> None:
    # todo
    switch_sol(sol_cb.get())

sol_cb.bind('<<ComboboxSelected>>', sol_changed)


switch_prob(int(prob_cb.get()))
switch_sol(sol_cb.get())


@key('Right')
def next_prob() -> None:
    name = int(prob_cb.get())
    pi = prob_ids.index(name)
    next_pi = (pi+1) % len(prob_ids)
    next_id = prob_ids[next_pi]
    prob_cb.set(next_id)
    switch_prob(next_id)


@key('Left')
def prev_prob() -> None:
    name = int(prob_cb.get())
    pi = prob_ids.index(name)
    next_pi = (pi-1) % len(prob_ids)
    next_id = prob_ids[next_pi]
    prob_cb.set(next_id)
    switch_prob(next_id)


@key('Next')
def next_sol() -> None:
    if cur_sol_ids == []:
        return
    name = sol_cb.get()
    pi = cur_sol_ids.index(name)
    next_pi = (pi+1) % len(cur_sol_ids)
    next_id = cur_sol_ids[next_pi]
    sol_cb.set(next_id)
    switch_sol(next_id)


@key('Prior')
def prev_sol() -> None:
    if cur_sol_ids == []:
        return
    name = sol_cb.get()
    pi = cur_sol_ids.index(name)
    next_pi = (pi-1) % len(cur_sol_ids)
    next_id = cur_sol_ids[next_pi]
    sol_cb.set(next_id)
    switch_sol(next_id)


@key('Escape')
def quit():
    exit(0)


legend_items = [f'{key}: {f.__name__}' for key, f in shortcuts.items()]
legend = ', '.join(legend_items)
bottom_label = tk.Label(root, text=(legend))
bottom_label.pack()


tk.mainloop()
