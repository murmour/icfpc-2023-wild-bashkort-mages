
import tkinter as tk
from tkinter import ttk
import sys
import 手


start_prob = sys.argv[1] if len(sys.argv) >= 2 else None
pixel_size = 1
suki_mode = False


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
sol_tags = []

sol_var = tk.StringVar()
sol_cb = ttk.Combobox(root, textvariable=sol_var)
sol_cb['values'] = []
sol_cb['state'] = 'readonly'
sol_cb.pack()


top_label = tk.Label(root)
top_label.pack()


canvas = tk.Canvas(root, width=0, height=0, cursor="crosshair")
canvas.pack(fill=tk.BOTH, expand=tk.YES)


coord_label = tk.Label(root)
coord_label.place(x=0, y=0)


# score_label = tk.Label(root)
# score_label.place(x=0, y=25)


# def set_score(score) -> None:
#     if score is None:
#         score_label.config(text=f'score: -')
#     else:
#         score_label.config(text=f'score: {score}')


def tk_color_from_rgba(rgba) -> str:
    (r, g, b, _) = rgba
    return f'#{r:02x}{g:02x}{b:02x}'


def project(x: float, y: float) -> tuple[float, float]:
    ch = canvas.winfo_height()
    xx = x*pixel_size
    yy = ch - (y*pixel_size)
    return (xx, yy)


def unproject(x: float, y: float) -> tuple[float, float]:
    # todo: handl_boxe zoom
    ch = canvas.winfo_height()
    y = ch-y
    return (x/pixel_size, y/pixel_size)


def on_mouse_move(event) -> None:
    (x, y) = unproject(event.x, event.y)
    coord_label.config(text=f'x: {x}, y: {y}')

canvas.bind("<Motion>", on_mouse_move)


def draw_dot(x, y, fill, tag) -> None:
    canvas.create_rectangle(
        x-1, y-1, x+1, y+1,
        fill=fill, width=0, tags=tag
    )


def fit_pixel_size(w: float, h: float) -> None:
    global pixel_size
    cw = canvas.winfo_width()
    ch = canvas.winfo_height()
    pixel_size = min(cw/w, ch/h)


def fit_pixel_size_to_prob() -> None:
    if prob is None:
        return
    rw = prob['room_width']
    rh = prob['room_height']
    fit_pixel_size(rw, rh)


def draw_prob() -> None:
    canvas.delete('prob')
    if suki_mode or prob is None:
        return
    fit_pixel_size_to_prob()

    # room
    rw = prob['room_width']
    rh = prob['room_height']
    (x0, y0) = project(0, 0)
    (x1, y1) = project(rw, rh)
    canvas.create_rectangle(
        x0, y0, x1, y1,
        fill='white', width=0, tags='prob'
    )

    # stage
    (sx, sy) = prob['stage_bottom_left']
    (x0, y0) = project(sx, sy)
    sw = prob['stage_width']
    sh = prob['stage_height']
    (x1, y1) = project(sx+sw, sy+sh)
    canvas.create_rectangle(
        x0, y0, x1, y1,
        fill='white', outline='red', width=1, tags='prob'
    )

    # pillars
    if 'pillars' in prob:
        for p in prob['pillars']:
            [ cx, cy ] = p['center']
            (cx_px, cy_px) = project(cx, cy)
            (r, _) = project(p['radius'], 0)
            canvas.create_oval(
                cx_px-r, cy_px-r, cx_px+r, cy_px+r,
                outline='blue', tags='prob'
            )

    # attendees
    for a in prob['attendees']:
        x = a['x']
        y = a['y']
        (x0, y0) = project(x, y)
        draw_dot(x0, y0, 'black', 'prob')


def draw_sol() -> None:
    canvas.delete('sol')
    if suki_mode or sol is None:
        return
    fit_pixel_size_to_prob()
    for p, v in zip(sol['placements'], sol['volumes']):
        x = p['x']
        y = p['y']
        (x0, y0) = project(x, y)
        (r, _) = project(5, 0)
        if v == 10:
            canvas.create_oval(x0-r, y0-r, x0+r, y0+r, tags='sol')
            continue
        if v == 0:
            canvas.create_oval(x0-r, y0-r, x0+r, y0+r, fill='grey', tags='sol')
            continue
        assert(False)


def suki_to_color(suki: float) -> str:
    if suki > 0:
        return 'red'
    if suki < 0:
        return 'blue'
    return 'white'


def draw_suki() -> None:
    canvas.delete('suki')
    if not suki_mode or prob is None:
        return

    ms = prob["musicians"]
    ats = prob["attendees"]
    ats_ct = len(ats)
    ins_ct = len(set(ms))
    fit_pixel_size(ats_ct, ins_ct)

    (x1, y1) = project(ats_ct, ins_ct)
    canvas.create_rectangle(
        0, 0, x1, y1,
        fill='white', width=0, tags='suki'
    )

    for j, a in enumerate(ats):
        for i in range(ins_ct):
            (x0, y0) = project(j, i)
            (x1, y1) = project(j+1, i+1)
            canvas.create_rectangle(
                x0, y0, x1, y1,
                fill = suki_to_color(a['tastes'][i]), width=0, tags='suki'
            )


def draw_all() -> None:
    draw_prob()
    draw_sol()
    draw_suki()


def on_resize(event) -> None:
    draw_all()

canvas.bind('<Configure>', on_resize)


def on_click1(event) -> None:
    xx = canvas.canvasx(event.x)
    yy = canvas.canvasy(event.y)
    canvas.scale('all', xx, yy, 1.1, 1.1)

canvas.bind("<Button-1>", on_click1)


def on_click2(event) -> None:
    draw_all()

canvas.bind("<Button-3>", on_click2)


def on_wheel_down(event):
    xx = canvas.canvasx(event.x)
    yy = canvas.canvasy(event.y)
    canvas.scale('all', xx, yy, 1.1, 1.1)


def on_wheel_up(event):
    xx = canvas.canvasx(event.x)
    yy = canvas.canvasy(event.y)
    factor = 1/1.1
    canvas.scale('all', xx, yy, factor, factor)


# root.bind("<MouseWheel>", mouse_wheel)
root.bind("<Button-5>", on_wheel_down)
root.bind("<Button-4>", on_wheel_up)


def update_top_label() -> None:
    p = prob
    score = f'{sol["score"]:,}' if sol is not None else '-'
    text = (
        f'musicians: {len(p["musicians"])}, '
        f'instuments: {len(set(p["musicians"]))}, '
        f'attendees: {len(p["attendees"])}, '
        f'pillars: {len(p["pillars"])}, '
        f'room_size: {p["room_width"]}x{p["room_height"]}, '
        f'stage_size: {p["stage_width"]}x{p["stage_height"]}, '
        f'score: {score}'
    )
    top_label.config(text=text)


def switch_sol(sol_tag: str) -> None:
    global sol
    if sol_tag == '':
        canvas.delete('sol')
        sol = None
        return
    if ':' in sol_tag:
        sol_tag = sol_tag.partition(':')[0]

    prob_id = int(prob_cb.get())
    sol = 手.get_sol(prob_id, sol_tag)
    draw_all()

    update_top_label()
    coord_label.focus_set()


def switch_prob(prob_id: int) -> None:
    global sol_tags, prob

    prob = 手.get_prob(prob_id)
    draw_all()

    best_score = 手.get_best_score(prob_id)
    def f(tag):
        return (tag, 手.get_score(prob_id, tag))
    sol_tags = [ f(tag) for tag in 手.get_sol_tags(prob_id) ]
    sol_tags.sort(key = lambda x: x[1], reverse = True)
    def f(score):
        if score == 0 or best_score == 0:
            return -1
        return score/best_score
    sol_tags = [ f'{x[0]}: {f(x[1]):0,.2f}' for x in sol_tags ]

    sol_cb['values'] = sol_tags
    if len(sol_tags) == 0:
        sol_cb.set('')
    else:
        sol_cb.current(0)
    switch_sol(sol_cb.get())

    update_top_label()
    coord_label.focus_set()


def prob_changed(event) -> None:
    switch_prob(int(prob_cb.get()))

prob_cb.bind('<<ComboboxSelected>>', prob_changed)


def sol_changed(event) -> None:
    switch_sol(sol_cb.get())

sol_cb.bind('<<ComboboxSelected>>', sol_changed)


switch_prob(int(prob_cb.get()))
switch_sol(sol_cb.get())


@key('s')
def toggle_suki() -> None:
    global suki_mode
    suki_mode = not suki_mode
    draw_all()


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
    if sol_tags == []:
        return
    name = sol_cb.get()
    pi = sol_tags.index(name)
    next_pi = (pi+1) % len(sol_tags)
    next_id = sol_tags[next_pi]
    sol_cb.set(next_id)
    switch_sol(next_id)


@key('Prior')
def prev_sol() -> None:
    if sol_tags == []:
        return
    name = sol_cb.get()
    pi = sol_tags.index(name)
    next_pi = (pi-1) % len(sol_tags)
    next_id = sol_tags[next_pi]
    sol_cb.set(next_id)
    switch_sol(next_id)


@key('Escape')
def quit() -> None:
    exit(0)


legend_items = [f'{key}: {f.__name__}' for key, f in shortcuts.items()]
legend = ', '.join(legend_items)
bottom_label = tk.Label(root, text=(legend))
bottom_label.pack()


tk.mainloop()
