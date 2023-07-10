
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

prob_count = 90


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


def get_prob(prob_id: int) -> dict:
    fname = get_prob_path(prob_id)
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


def save_prob(prob_id: int, prob: dict) -> None:
    fname = get_prob_path(prob_id)
    with io.open(fname, 'w') as f:
        f.write(json.dumps(prob))


def download_and_save_prob(prob_id: int) -> None:
    print(f'Downloading problem {prob_id}...')
    prob = get_prob_request(prob_id)
    save_prob(prob_id, prob)
    print(f'Saved {prob_id}')
    sol_dir = get_sol_dir_path(prob_id)
    if not os.path.exists(sol_dir):
        os.mkdir(sol_dir)


# Make a submission. Returns a submission ID.
def send_sol_request(prob_id: int, sol_tag: str) -> str:
    print(f'Sending {prob_id}.{sol_tag}...')
    sol = get_sol(prob_id, sol_tag)
    if 'score' in sol:
        del sol['score']
    prob = get_prob(prob_id)
    nmus = len(prob['musicians'])
    if 'volumes' not in sol:
        print(f'{prob_id}.{sol_tag} has no volumes!')
        sol['volumes'] = [10.0] * nmus

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


def get_sol(prob_id: int, tag: str) -> dict | None:
    sol_dir = get_sol_dir_path(prob_id)
    fname = f'{sol_dir}/{tag}.solution'
    if not os.path.exists(fname):
        return None
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


def get_score(prob_id: int, sol_tag: str) -> float:
    sol = get_sol(prob_id, sol_tag)
    score = sol['score']
    if 'volumes' not in sol:
        print(f'{prob_id}.{sol_tag} has no volumes!')
        score *= 10
    return score


def score_sol_by_tag(prob_id: int, sol_tag: str) -> float:
    prob_path = get_prob_path(prob_id)
    sol_dir = get_sol_dir_path(prob_id)
    sol_path = f'{sol_dir}/{sol_tag}.solution'
    sys.stdout.flush()
    try:
        sol = subprocess.check_output([
            '脳/脳',
            '-pp', prob_path,
            '-score', sol_path,
        ])
        return float(sol.decode())
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print(ex)
        sys.stdout.flush()
        return 0


def score_sol(prob_id: int, sol: dict) -> float:
    prob_path = get_prob_path(prob_id)
    sys.stdout.flush()
    try:
        sol = subprocess.check_output([
                '脳/脳',
                '-pp', prob_path,
                '-score', 'stdin',
            ],
            input=json.dumps(sol).encode()
        )
        return float(score.decode())
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print(ex)
        sys.stdout.flush()
        return 0


def get_best_sol_tag(prob_id) -> str | None:
    best_score = 0
    best_tag = None
    for tag in get_sol_tags(prob_id):
        score = get_score(prob_id, tag)
        if score > best_score:
            best_score = score
            best_tag = tag
    print(f'{prob_id}: {best_tag} ({best_score:,})')
    return best_tag


def send_all_best() -> dict:
    for prob_id in all_prob_ids():
        sol_tag = get_best_sol_tag(prob_id)
        if sol_tag is None:
            print(f'{prob_id}: -')
            continue
        send_sol_request(prob_id, sol_tag)


def send_all_best_lazy(userboard) -> dict:
    for prob_id in all_prob_ids():
        sol_tag = get_best_sol_tag(prob_id)
        if sol_tag is None:
            print(f'{prob_id}: -')
            continue
        score = get_score(prob_id, sol_tag)
        remote_score = userboard["problems"][prob_id-1]
        if remote_score is None or score > remote_score:
            print(f'{prob_id}: {score:,} > {remote_score:,}')
            send_sol_request(prob_id, sol_tag)


def print_prob_stats(prob_id: int) -> dict:
    p = get_prob(prob_id)
    print(f'{prob_id}:')
    print(f'  musicians: {len(p["musicians"])}')
    print(f'  instruments: {len(set(p["musicians"]))}')
    print(f'  attendees: {len(p["attendees"])}')
    print(f'  pillars: {len(p["pillars"])}')
    print(f'  room_size: {p["room_width"]}x{p["room_height"]}')
    print(f'  stage_size: {p["stage_width"]}x{p["stage_height"]}')


def print_sol_stats(print_all_sols = True) -> None:
    total_score = 0
    for prob_id in all_prob_ids():
        if print_all_sols:
            print(f'{prob_id}:')
        sol_tags = get_sol_tags(prob_id)
        best = None
        best_score = 0
        for tag in sol_tags:
            score = get_score(prob_id, tag)
            if print_all_sols:
                print(f'  {tag} = {score:,}')
            if best is None or score > best_score:
                best = tag
                best_score = score

        if best is not None:
            if print_all_sols:
                print(f'  BEST = {best}({best_score:,})')
            else:
                print(f'{prob_id}: {best}({best_score:,})')
            total_score += max(best_score, 0)
    print()
    print(f'TOTAL = {total_score:,}')


# solver CLI:
# -p 1  // problem id
# -pp path  // problem path
# -out filename // output file name (default=stdout)
# -timeout 120 // timeout in seconds (default = 120)
# -s two_row // solver (two_row/border/regular)
# -score solution_file_name // score solution and print the score to stdout (doesn't run the solver)


def get_solver_cmd(solver_id: str, solver_args: str, prob_id: int) -> list[str]:
    prob_path = get_prob_path(prob_id)

    # problem has mask?
    prob = get_prob(prob_id)
    mask = None
    if 'mask' in prob:
        mask = prob['mask']

    print(f'{solver_id}({prob_id})...')
    cmd = ['脳/脳']
    cmd += ['-s', solver_id]
    cmd += ['-pp', prob_path]
    if solver_args != '':
        cmd += solver_args.split(' ')
    if mask is not None:
        cmd += ['-mask', str(mask)]

    return cmd


def run_solver(solver_id: str, solver_args: str, prob_id: int) -> str | None:
    cmd = get_solver_cmd(solver_id, solver_args, prob_id)
    sys.stdout.flush()
    try:
        return subprocess.check_output(cmd).decode()
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print('OK')
        print(ex)
        sys.stdout.flush()
        return None


def save_sol(solver_id: str, prob_id: int, sol_tag: str, sol: dict) -> None:
    print(f'saving {solver_id}({prob_id}) as {sol_tag}')
    sol_dir = get_sol_dir_path(prob_id)
    sol_path = f'{sol_dir}/{sol_tag}.solution'
    with io.open(sol_path, 'wb') as h:
        h.write(json.dumps(sol).encode())


def solve(solver_id: str, solver_args: str, sol_tag: str, prob_id: int) -> None:
    sol_str = run_solver(solver_id, solver_args, prob_id)
    if sol_str is None:
        print(f'{solver_id}({prob_id}) -> ERROR')
        return
    sol = json.loads(sol_str)
    score = sol['score']
    print(f'{solver_id}({prob_id}) -> OK ({score:,})')
    save_sol(solver_id, prob_id, sol_tag, sol)


def parsolve(
        solver_id: str,
        solver_args: str,
        sol_tag: str,
        prob_ids: list,
        job_count: int
) -> None:

    def start_solving(prob_id: int) -> dict:
        cmd = get_solver_cmd(solver_id, solver_args, prob_id)
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE
        )
        return { 'prob_id': prob_id, 'process': process }

    def finish_solving(job, retcode: int) -> None:
        prob_id = job['prob_id']
        if retcode != 0:
            print(f'{solver_id}({prob_id}) -> ERROR ({retcode})')
            return
        sol_str = job['process'].communicate()[0].decode()
        sol = json.loads(sol_str)
        score = sol['score']
        print(f'{solver_id}({prob_id}) -> OK ({score:,})')
        save_sol(solver_id, prob_id, sol_tag, sol)

    queue = list(prob_ids)
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
    sol_tag = get_best_sol_tag(prob_id)
    return 0 if sol_tag is None else get_score(prob_id, sol_tag)


def get_sorted_scores():
    data = []
    for pid in all_prob_ids():
        data.append((get_best_score(pid), pid))
    data.sort(reverse=True)
    tot = 0
    for score, pid in data:
        tot += max(score, 0)
        print('%12.0f %3d' % (score, pid))
    print('Total = %.0f', tot)


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
        print(f'{i+1}. {line["username"]} = {line["score"]:,}')


def print_userboard(ub: dict) -> None:
    print(f'userboard:')
    for i, score in enumerate(ub['problems']):
        print(f'{i+1} = {(0 if score is None else score):,}')


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


def purge_bad_sols() -> None:
    for prob_id in all_prob_ids():
        sol_dir = get_sol_dir_path(prob_id)
        for tag in get_sol_tags(prob_id):
            score = get_score(prob_id, tag)
            if score == 0:
                print(f'bad: {prob_id}.{tag}')
                sol_path = f'{sol_dir}/{tag}.solution'
                # эту строку нужно раскомментировать при каждом использовании:
                # os.remove(sol_path)


def patch_scores() -> None:
    for prob_id in all_prob_ids():
        for sol_tag in get_sol_tags(prob_id):
            sol = get_sol(prob_id, sol_tag)
            if 'score' not in sol:
                score = score_sol_by_tag(prob_id, sol_tag)
                if score == 0:
                    print(f'bad sol: {prob_id}.{sol_tag}')
                    continue
                print(f'adding score: {prob_id}.{sol_tag} = {score:,}')
                sol['score'] = score
                # эту строку нужно раскомментировать при каждом использовании:
                # save_sol('patch_scores', prob_id, sol_tag, sol)


def patch_problems_with_masks() -> None:
    for prob_id in all_prob_ids():
        sol = get_sol(prob_id, 'smart_1000_probe_c')
        if sol is None:
            print(f'{prob_id}: -')
            continue
        mask = sol['mask']
        prob = get_prob(prob_id)
        prob['mask'] = mask
        print(f'{prob_id}: {mask}')
        if mask != -1:
            # эту строку нужно раскомментировать при каждом использовании:
            save_prob(prob_id, prob)


def volumize_sol(prob_id: int, sol: dict) -> dict:
    prob_path = get_prob_path(prob_id)
    sys.stdout.flush()
    try:
        sol = subprocess.check_output([
                '脳/脳',
                '-pp', prob_path,
                '-volumize', 'stdin',
            ],
            input=json.dumps(sol).encode()
        )
        return json.loads(sol.decode())
    except subprocess.CalledProcessError as ex: # error code <> 0
        print('-------- ERROR --------')
        print(ex)
        sys.stdout.flush()
        return None


def volumize_all_sols() -> None:
    for prob_id in all_prob_ids():
        for sol_tag in get_sol_tags(prob_id):
            sol = get_sol(prob_id, sol_tag)
            score = sol['score']
            vol_sol = volumize_sol(prob_id, sol)
            if vol_sol is None:
                continue
            vol_score = vol_sol['score']
            if vol_score == score:
                print(f'{prob_id}.{sol_tag}: same')
            elif vol_score < score:
                print(f'{prob_id}.{sol_tag}: lower!')
            elif vol_score > score:
                win = vol_score - score*10
                print(f'{prob_id}.{sol_tag}: ratio: {vol_score/score}, win: {win}')
                # эту строку нужно раскомментировать при каждом использовании:
                # save_sol('volumize', prob_id, sol_tag, vol_sol)


def dp_is_best(prob_id: int) -> bool:
    sol_tag = get_best_sol_tag(prob_id)
    if sol_tag is None:
        return False
    for pre in ['dp_', 'dpf_', 'dpff_', 'dpass_']:
        if sol_tag.startswith(pre):
            return True
    print(f'{prob_id}\'s best is not dp')
    return False


if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'solve':
        solver = sys.argv[2]
        solver_args = sys.argv[3]
        sol_id = sys.argv[4]
        start_id = int(sys.argv[5])
        end_id = int(sys.argv[6])
        for prob_id in range(start_id, end_id+1):
            solve(solver, solver_args, sol_id, prob_id)
        exit(0)

    if cmd == 'parsolve':
        solver_id = sys.argv[2]
        solver_args = sys.argv[3]
        sol_id = sys.argv[4]
        start_id = int(sys.argv[5])
        end_id = int(sys.argv[6])
        job_count = int(sys.argv[7])
        prob_ids = range(start_id, end_id+1)
        if solver_id == 'dp_if_best':
            solver_id = 'dp'
            prob_ids = [ pid for pid in prob_ids if dp_is_best(pid) ]
        parsolve(solver_id, solver_args, sol_id, prob_ids, job_count)
        exit(0)

    if cmd == 'download_probs':
        for pid in all_prob_ids():
            download_and_save_prob(pid)
        exit(0)

    if cmd == 'get_prob_stats':
        for pid in all_prob_ids():
            print_prob_stats(pid)
            print()
        exit(0)

    if cmd == 'get_sol_stats':
        print_sol_stats(print_all_sols = True)
        exit(0)

    if cmd == 'get_best_sols':
        print_sol_stats(print_all_sols = False)
        exit(0)

    if cmd == 'get_sorted_scores':
        get_sorted_scores()
        exit(0)

    if cmd == 'send_all_best':
        send_all_best()
        exit(0)

    if cmd == 'send_all_best_lazy':
        ub = get_userboard_request()
        send_all_best_lazy(ub)
        exit(0)

    if cmd == 'purge_bad_sols':
        purge_bad_sols()
        exit(0)

    if cmd == 'patch_scores':
        patch_scores()
        exit(0)

    if cmd == 'patch_problems_with_masks':
        patch_problems_with_masks()
        exit(0)

    if cmd == 'volumize_sols':
        volumize_all_sols()
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
