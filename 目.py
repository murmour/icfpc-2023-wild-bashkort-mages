
import tkinter as tk
from tkinter import ttk
import sys
import 手


start_prob = sys.argv[1] if len(sys.argv) >= 2 else None
pixel_size = 1


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


prob = None
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


sol = None
sol_ids = []

sol_var = tk.StringVar()
sol_cb = ttk.Combobox(root, textvariable=sol_var)
sol_cb['values'] = []
sol_cb['state'] = 'readonly'
sol_cb.pack()


top_label = tk.Label(root)
top_label.pack()


canvas = tk.Canvas(root, width=0, height=0)
canvas.pack(fill=tk.BOTH, expand=tk.YES)


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


def project_coord(x: float, y: float) -> tuple[float, float]:
    ch = canvas.winfo_height()
    return (x*pixel_size, ch - (y*pixel_size))


def on_mouse_move(event) -> None:
    (x, y) = project_coord(event.x, event.y)
    coord_label.config(text=f'x: {x}, y: {y}')

canvas.bind("<Motion>", on_mouse_move)


def draw_prob(prob) -> None:
    canvas.delete('prob')
    if prob is None:
        return

    # room
    rw = prob['room_width']
    rh = prob['room_height']
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
        canvas.create_rectangle(
            x0, y0, x0+2, y0+2,
            fill='black', width=0, tags='prob'
        )


def draw_sol(sol) -> None:
    canvas.delete('sol')
    if sol is None:
        return
    for m in sol['placements']:
        x = m['x']
        y = m['y']
        (x0, y0) = project_coord(x, y)
        (x1, y1) = project_coord(x+1, y+1)
        canvas.create_rectangle(
            x0, y0, x1, y1,
            fill='black', width=0, tags='sol'
        )


def fit_pixel_size(w: float, h: float) -> None:
    global pixel_size
    cw = canvas.winfo_width()
    ch = canvas.winfo_height()
    pixel_size = min(cw/w, ch/h)


def fit_pixel_size_to_prob(prob) -> None:
    if prob is None:
        return
    rw = prob['room_width']
    rh = prob['room_height']
    fit_pixel_size(rw, rh)


def on_resize(event):
    old_pixel_size = pixel_size
    fit_pixel_size_to_prob(prob)
    if pixel_size != old_pixel_size:
        draw_prob(prob)
        draw_sol(sol)

canvas.bind('<Configure>', on_resize)


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


def switch_sol(sol_tag: str) -> None:
    global sol
    if sol_tag == '':
        return
    prob_id = int(prob_cb.get())
    sol = 手.get_sol(prob_id, sol_tag)
    draw_sol(sol)

    coord_label.focus_set()


def switch_prob(prob_id: int) -> None:
    global sol_ids, prob

    prob = 手.get_prob(prob_id)
    fit_pixel_size_to_prob(prob)
    draw_prob(prob)
    update_top_label(prob)

    sol_ids = 手.get_sol_tags(prob_id)
    sol_cb['values'] = sol_ids
    if len(sol_ids) == 0:
        sol_cb.set('')
    else:
        sol_cb.current(0)
    switch_sol(sol_cb.get())

    coord_label.focus_set()


def prob_changed(event) -> None:
    switch_prob(int(prob_cb.get()))

prob_cb.bind('<<ComboboxSelected>>', prob_changed)


def sol_changed(event) -> None:
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
    if sol_ids == []:
        return
    name = sol_cb.get()
    pi = sol_ids.index(name)
    next_pi = (pi+1) % len(sol_ids)
    next_id = sol_ids[next_pi]
    sol_cb.set(next_id)
    switch_sol(next_id)


@key('Prior')
def prev_sol() -> None:
    if sol_ids == []:
        return
    name = sol_cb.get()
    pi = sol_ids.index(name)
    next_pi = (pi-1) % len(sol_ids)
    next_id = sol_ids[next_pi]
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
