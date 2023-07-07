
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

prob_count = 45


def all_prob_ids() -> list[int]:
    return [ i+1 for i in range(prob_count) ]


def get_prob_request(id: int) -> dict:
    res = subprocess.check_output(
        ['curl',
         '--compressed',
         '-H', f'Authorization: Bearer {api_key}',
         f'{api}/problem?problem_id={id}' ])
    time.sleep(1)
    res_js = json.loads(res.decode('utf-8'))
    assert('Success' in res_js)
    return json.loads(res_js['Success'])


def get_prob_path(id: int) -> str:
    path = f'問/{id}.problem'


def get_sol_dir_path(prob_id: int) -> str:
    path = f'答/{prob_id}'
    assert(os.path.exists(path))
    return path


def download_and_save_problem(id: int):
    print(f'Downloading problem {id}...')
    js = get_prob_request(id)
    fname = get_prob_path(id)
    with io.open(fname, 'w') as f:
        f.write(json.dumps(js))
        print(f'Saved to {fname}')


def get_prob(id: int) -> dict:
    fname = get_prob_path(id)
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


def get_sol_tags(id: int) -> list:
    sols_dir = get_sol_dir_path(id)
    tags = [ f[:-9] for f in os.listdir(sols_dir) if f.endswith('.solution') ]
    tags.sort()
    return tags


def get_sol(id: int, tag: str) -> dict:
    sols_dir = get_sol_dir_path(id)
    fname = f'{sols_dir}/{tag}.solution'
    with io.open(fname, 'r') as h:
        return json.loads(h.read())


def print_prob_stats(id: int) -> dict:
    p = get_prob(id)
    print(f'{id}.problem:')
    print(f'musicians: {len(p["musicians"])}')
    print(f'instuments: {len(set(p["musicians"]))}')
    print(f'attendees: {len(p["attendees"])}')
    print(f'room_size: {p["room_width"]}x{p["room_height"]}')
    print(f'stage_size: {p["stage_width"]}x{p["stage_height"]}')
    print()


if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'get_probs':
        for pid in all_prob_ids():
            download_and_save_problem(pid)
        exit(0)

    if cmd == 'print_prob_stats':
        for pid in all_prob_ids():
            print_prob_stats(pid)
        exit(0)

    print(f'invalid command: {cmd}')
    exit(1)
