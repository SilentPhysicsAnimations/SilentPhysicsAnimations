import os
import shutil
import subprocess
import time

# Ask for desired animations

animations = [
    { 'name': "Projectile Motion",        'file': "01-ProjectileMotion.py"        },
    { 'name': "Inclined Plane",           'file': "02-InclinedPlane.py"           },
    { 'name': "Centripetal Acceleration", 'file': "03-CentripetalAcceleration.py" },
    { 'name': "Planetary Motion",         'file': "04-PlanetaryMotion.py"         },
    { 'name': "Ballistic Pendulum",       'file': "05-BallisticPendulum.py"       },
    { 'name': "Archimedes' Principle",    'file': "06-ArchimedesPrinciple.py"     },
    { 'name': "Simple Harmonic Motion",   'file': "07-SimpleHarmonicMotion.py"    },
    { 'name': "Doppler Effect",           'file': "08-DopplerEffect.py"           },
    { 'name': "Ideal Gas Law",            'file': "09-IdealGasLaw.py"             },
    { 'name': "Electric Field Lines",     'file': "10-ElectricFieldLines.py"      },
    { 'name': "Velocity Selector",        'file': "11-VelocitySelector.py"        },
    { 'name': "Time Dilation",            'file': "12-TimeDilation.py"            },
]

print("[33m⣿⣿⣿⣿ Which animations do you want to render?[0m")

for i, anim in enumerate(animations, start=1):
    name = anim['name']
    print(f"[1;2m[{i:02}][0m {name}", end="\n" if i % 2 == 0 else (26 - len(name))*" ")

while True:
    render_numbers = []
    raw_input = input("[33m⣿⣿⣿⣿ (input number, list, or range; ex. `9`, `1,2,5`, `1-12`) [0m")
    if raw_input == '':
        raw_input = '1-12'
    parts = raw_input.replace(' ', '').split(',')
    try:
        for part in parts:
            ends = list(map(int, part.split('-', 1)))
            if len(ends) == 1:
                ends.append(ends[0])
            if not (1 <= ends[0] <= 12 and 1 <= ends[1] <= 12):
                raise ValueError
            for num in range(ends[0], ends[1]+1):
                if num not in render_numbers:
                    render_numbers.append(num)
        break
    except ValueError:
        print("[31m⣿⣿⣿⣿ Sorry, I don't understand. Please try again.[0m")
        
animations_to_render = list(map(lambda n: animations[n-1], render_numbers))

# Ask for desired quality

print(
    "[33m⣿⣿⣿⣿ In what quality do you want to render?[0m",
    "[1;2m[01][0m HD/720p @ 30fps",
    "[1;2m[02][0m FHD/1080p @ 60fps [2m(recommended)[0m",
    "[1;2m[03][0m QHD/2K/1440p @ 60fps",
    "[1;2m[04][0m UHD/4K/2160p @ 60fps",
    sep='\n'
)

while True:
    raw_input = input("[33m⣿⣿⣿⣿ (input number) [0m")
    if raw_input == '':
        raw_input = '2'
    try:
        quality_number = int(raw_input)
        if not 1 <= quality_number <= 4:
            raise ValueError
        break
    except ValueError:
        print("[31m⣿⣿⣿⣿ Sorry, I don't understand. Please try again.[0m")

quality = ['m', 'h', 'p', 'k'][quality_number-1]
quality_folder = ['720p30', '1080p60', '1440p60', '2160p60'][quality_number-1]

# Start rendering

spinner = ['⣿⣿⣿⠀', '⢸⣿⣿⡇', '⠀⣿⣿⣿', '⠀⢸⣿⣿', '⠀⠀⣿⣿', '⠀⠀⢸⣿', '⠀⠀⠀⣿', '⠀⠀⠀⢸', '⠀⠀⠀⠀', '⠀⠀⠀⠀', '⡇⠀⠀⠀', '⣿⠀⠀⠀', '⣿⡇⠀⠀', '⣿⣿⠀⠀', '⣿⣿⡇⠀']

spinner_index = 0
spinner_max = len(spinner)

print('[?25l', end='')

for anim in animations_to_render:
    name = anim['name']
    print(f"[1;2m⠀⠀⠀⠀ [00:00:00][0;2m     Waiting to render [1;2m‘{name}’[0m")

num_anims = len(animations)

render_start = time.time()

for line, anim in enumerate(animations_to_render):
    name, file = anim['name'], anim['file']
    base = os.path.splitext(file)[0]

    file_render_start = time.time()

    render_process = subprocess.Popen([
        'manim', '-o', f'{base}.mp4', '--format', 'mp4', '--disable_caching', '-q', quality,
        file, f'{base[3:]}Scene'
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    while True:
        try:
            return_code = render_process.wait(0.1)
            now = time.time()
            elapsed = time.strftime('%X', time.gmtime(now - file_render_start))
            if return_code == 0:
                print(
                    f'[{num_anims - line}F[K',
                    f"[1;32m⣿⣿⣿⣿ [{elapsed}][0m Successfully rendered [1;32m‘{name}’[0m",
                    f'[{num_anims - line}E',
                    sep='', end='', flush=True
                )
                shutil.copy2(
                    os.path.join('media', 'videos', base, quality_folder, f'{base}.mp4'),
                    os.path.join('videos', quality_folder, f'{base}.mp4')
                )
            else:
                print(
                    f'[{num_anims - line}F[K',
                    f"[1;31m⣿⣿⣿⣿ [{elapsed}][0;31m      Could not render [1;31m‘{name}’[0m",
                    f'[{num_anims - line}E',
                    sep='', end='', flush=True
                )
            break
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            now = time.time()
            elapsed = time.strftime('%X', time.gmtime(now - file_render_start))
            print(
                f'[{num_anims - line}F[K',
                f"[1;34m{spinner[spinner_index]} [{elapsed}][0m   Rendering animation [1;34m‘{name}’[0m",
                f'[{num_anims - line}E',
                sep='', end='', flush=True
            )
            total_elapsed = time.strftime('%X', time.gmtime(now - render_start))
            print(
                '[G[K',
                f"[1;33m{spinner[spinner_index]} [{total_elapsed}][0m Rendering animations into [1;33m`videos/{quality_folder}/`[0m",
                sep='', end='', flush=True
            )
            spinner_index = (spinner_index + 1) % spinner_max
            continue

print(
    '[G[K',
    f"[1;33m⣿⣿⣿⣿ [{total_elapsed}][0m Finished rendering animations into [1;33m`videos/{quality_folder}/`[0m",
    sep='', end='', flush=True
)

print('[?25h')
