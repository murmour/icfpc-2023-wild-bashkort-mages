
import os
import io
import re
import sys
import subprocess
import json
import time
import shutil


api = 'https://api.icfpcontest.com'
api_key = (
    'eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJ1aWQiOiI2NDlkZTZhMjh'
    'jNjg1MzEzZDFjNjBkNjEiLCJpYXQiOjE2ODgwNjk3OTQsImV4cCI6MTY5ODA'
    '2OTc5NH0.svNmE8_E5X7tEQjnqNF9rtdlgOckN21pZJ85br6uFOw'
)

prob_count = 55


def all_prob_ids() -> list[int]:
    return [ i+1 for i in range(prob_count) ]


def get_prob_path(prob_id: int) -> str:
    return f'問/{prob_id}.problem'


def get_sol_dir_path(prob_id: int) -> str:
    return f'答/{prob_id}'


def get_prob_request(prob_id: int) -> dict:
    res = subprocess.check_output([
        'curl',
        '--compressed',
        '-H', f'Authorization: Bearer {api_key}',
        f'{api}/problem?problem_id={prob_id}'
    ])
    time.sleep(1)
    res_js = json.loads(res.decode('utf-8'))
    assert('Success' in res_js)
    return json.loads(res_js['Success'])


def download_and_save_problem(prob_id: int):
    print(f'Downloading problem {prob_id}...')
    js = get_prob_request(prob_id)
    fname = get_prob_path(prob_id)
    with io.open(fname, 'w') as f:
        f.write(json.dumps(js))
        print(f'Saved to {fname}')
    sol_dir = get_sol_dir_path(prob_id)
    if not os.path.exists(sol_dir):
        os.mkdir(sol_dir)


def get_prob(prob_id: int) -> dict:
    fname = get_prob_path(prob_id)
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


# Make a submission. Returns a submission ID.
def send_sol_request(prob_id: int, sol_tag: str) -> str:
    print(f'Sending {prob_id}.{sol_tag}...')
    sol = get_sol(prob_id, sol_tag)
    sub = { 'problem_id': prob_id, 'contents': json.dumps(sol) }
    sub_path = 'submission.js'
    with io.open(sub_path, 'w') as f:
        f.write(json.dumps(sub))
    res = subprocess.check_output([
        'curl',
        '--json', f'@{sub_path}',
        '-H', f'Authorization: Bearer {api_key}',
        f'{api}/submission',
    ])
    time.sleep(1)
    res = res.decode('utf-8')
    print(f'response: {res}')
    return res


def get_sol_tags(prob_id: int) -> list:
    sol_dir = get_sol_dir_path(prob_id)
    tags = [ f[:-9] for f in os.listdir(sol_dir) if f.endswith('.solution') ]
    tags.sort()
    return tags


def get_sol(prob_id: int, tag: str) -> dict:
    sol_dir = get_sol_dir_path(prob_id)
    fname = f'{sol_dir}/{tag}.solution'
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


def get_score(prob_id: int, sol_tag: str) -> float | None:
    prob_path = get_prob_path(prob_id)
    sol_dir = get_sol_dir_path(prob_id)
    sol_path = f'{sol_dir}/{sol_tag}.solution'
    sys.stdout.flush()
    try:
        sol = subprocess.check_output([
            '../solver/solver',
            '-pp', prob_path,
            '-score', sol_path,
        ])
        return float(sol.decode())
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print(ex)
        sys.stdout.flush()
        return None


def get_best_sol(prob_id) -> str | None:
    best_score = None
    best_tag = None
    for tag in get_sol_tags(prob_id):
        score = get_score(prob_id, tag)
        if score is None:
            continue
        if best_score == None or score > best_score:
            best_score = score
            best_tag = tag
    if best_score == None:
        print(f'{prob_id}: -')
        return None
    print(f'{prob_id}: {best_tag} ({score})')
    return best_tag


def get_local_best() -> dict:
    r = {}
    for prob_id in all_prob_ids():
        sol = get_best_sol(prob_id)
        if sol is not None:
            r[prob_id] = sol
    return r


def print_prob_stats(prob_id: int) -> dict:
    p = get_prob(prob_id)
    print(f'{prob_id}:')
    print(f'  musicians: {len(p["musicians"])}')
    print(f'  instuments: {len(set(p["musicians"]))}')
    print(f'  attendees: {len(p["attendees"])}')
    print(f'  room_size: {p["room_width"]}x{p["room_height"]}')
    print(f'  stage_size: {p["stage_width"]}x{p["stage_height"]}')


total_score = 0

def print_sol_stats(prob_id: int) -> dict:
    global total_score

    print(f'{prob_id}:')
    sol_tags = get_sol_tags(prob_id)

    best = None
    best_score = None
    for tag in sol_tags:
        score = get_score(prob_id, tag)
        print(f'  {tag} = {score}')
        if best is None or score > best_score:
            best = tag
            best_score = score

    if best is not None:
        print(f'  BEST = {best}({best_score})')
        total_score += best_score


# solver CLI:
# -p 1  // problem id
# -pp path  // problem path
# -timeout 120 // timeout in seconds (default = 120)
# -out filename // output file name (default=stdout)
# -solver two_row // solver (two_row/border/regular)
# -score solution_file_name // score solution and print the score to stdout (doesn't run the solver)


def get_solver_cmd(solver_id: str, solver_args: str, prob_id: int) -> list[str]:
    prob_path = get_prob_path(prob_id)
    print(f'{solver_id}({prob_id})...')
    cmd = ['../solver/solver']
    cmd += ['-solver', solver_id]
    cmd += ['-pp', prob_path]
    if solver_args != '':
        cmd += solver_args.split(' ')
    return cmd


def run_solver(solver_id: str, solver_args: str, prob_id: int) -> str | None:
    cmd = get_solver_cmd(solver_id, solver_args, prob_id)
    sys.stdout.flush()
    try:
        sol = subprocess.check_output(cmd)
        return sol.decode()
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print('OK')
        print(ex)
        sys.stdout.flush()
        return None


def save_sol(solver_id: str, prob_id: int, sol_tag: str, sol: str) -> None:
    print(f'saving {solver_id}({prob_id}) as {sol_tag}')
    sol_dir = get_sol_dir_path(prob_id)
    sol_path = f'{sol_dir}/{sol_tag}.solution'
    with io.open(sol_path, 'wb') as h:
        h.write(sol.encode())


def solve(solver_id: str, solver_args: str, sol_tag: str, prob_id: int) -> None:
    sol = run_solver(solver_id, solver_args, prob_id)
    if sol is None:
        print(f'{solver_id}({prob_id}) -> ERROR')
        return
    sol_js = json.loads(sol)
    score = get_score(sol_js)
    print(f'{solver_id}({prob_id}) -> OK ({score})')
    save_sol(solver_id, prob_id, sol_tag, sol)


def parsolve(solver_id, solver_args, sol_tag, start_id, end_id, job_count) -> None:

    def start_solving(prob_id: int) -> dict:
        cmd = get_solver_cmd(solver_id, solver_args, prob_id)
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE
        )
        return { 'prob_id': prob_id, 'process': process }

    def finish_solving(job, retcode) -> None:
        prob_id = job['prob_id']
        if retcode != 0:
            print(f'{solver_id}({prob_id}) -> ERROR ({retcode})')
            return
        sol = job['process'].communicate()[0].decode()
        score = get_score(sol)
        print(f'{solver_id}({prob_id}) -> OK ({score})')
        save_sol(solver_id, prob_id, sol_tag, sol)

    queue = list(range(start_id, end_id+1))
    pool = [None] * job_count
    left = len(queue)
    while left > 0:
        time.sleep(0.01)
        for i in range(len(pool)):
            if pool[i] is None:
                if queue == []:
                    continue
                p = queue[0]
                queue = queue[1:]
                pool[i] = start_solving(p)
                if pool[i] is None:
                    left -= 1
                    continue
                else:
                    break
            job = pool[i]
            retcode = job['process'].poll()
            if retcode is not None:
                left -= 1
                pool[i] = None
                finish_solving(job, retcode)


def get_best_score(prob_id: int) -> float:
    sol_tags = get_sol_tags(prob_id)
    best = None
    best_score = None
    for tag in sol_tags:
        score = get_score(prob_id, tag)
        if best is None or score > best_score:
            best = tag
            best_score = score
    return best_score


def get_sorted_scores():
    data = []
    for pid in all_prob_ids():
        data.append((get_best_score(pid), pid))
    data.sort(reverse=True)
    for score, pid in data:
        print('%12.0f %3d' % (score, pid))


def get_scoreboard_request() -> dict:
    res = subprocess.check_output([
        'curl',
        '--compressed',
        f'{api}/scoreboard'
    ])
    time.sleep(1)
    return json.loads(res.decode('utf-8'))


def get_userboard_request() -> dict:
    res = subprocess.check_output([
        'curl',
        '--compressed',
        '-H', f'Authorization: Bearer {api_key}',
        f'{api}/userboard'
    ])
    time.sleep(1)
    res_js = json.loads(res.decode('utf-8'))
    assert('Success' in res_js)
    return res_js['Success']


def print_scoreboard(sb: dict) -> None:
    print(f'scoreboard at {sb["updated_at"]}:')
    for i, line in enumerate(sb['scoreboard']):
        print(f'{i}. {line["username"]} = {line["score"]}')


def print_userboard(sb: dict) -> None:
    print(f'userboard:')
    for i, score in enumerate(sb['problems']):
        print(f'{i} = {"-" if score is None else score}')


def update_username_request(username: str) -> None:
    data = json.dumps({ 'username': username })
    res = subprocess.check_output([
        'curl',
        '--compressed',
        '-H', f'Authorization: Bearer {api_key}',
        f'{api}/user/update_username',
        '--json', data
    ])
    time.sleep(1)
    print(res)
    res_js = json.loads(res.decode('utf-8'))
    assert('Success' in res_js)


if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'download_probs':
        for pid in all_prob_ids():
            download_and_save_problem(pid)
        exit(0)

    if cmd == 'get_prob_stats':
        for pid in all_prob_ids():
            print_prob_stats(pid)
            print()
        exit(0)

    if cmd == 'get_sol_stats':
        for pid in all_prob_ids():
            print_sol_stats(pid)
            print()
        print(f'TOTAL = {total_score}')
        exit(0)

    if cmd == 'get_sorted_scores':
        get_sorted_scores()
        exit(0)

    if cmd == 'send_all_best':
        best = get_local_best()
        for prob_id, sol_tag in best.items():
            send_sol_request(prob_id, sol_tag)
        exit(0)

    if cmd == 'get_scoreboard':
        sb = get_scoreboard_request()
        print_scoreboard(sb)
        exit(0)

    if cmd == 'get_userboard':
        ub = get_userboard_request()
        print_userboard(ub)
        exit(0)

    if cmd == 'update_username':
        update_username_request("WILD BASHKORT MAGES")
        exit(0)

    print(f'invalid command: {cmd}')
    exit(1)
