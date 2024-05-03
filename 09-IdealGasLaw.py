from manim import *

##########################

PRESSURE_COLOR = BLUE
VOLUME_COLOR = BLUE

NORMAL_TEMP_COLOR = YELLOW
LOW_TEMP_COLOR = BLUE_A
HIGH_TEMP_COLOR = RED

##########################

class Lock(VGroup):
    def __init__(self, shackle_width=0.2, shackle_radius=0.6, shackle_height=1, color=WHITE, **kwargs):
        body = RoundedRectangle(width=1.618, height=1, corner_radius=0.05, fill_opacity=1)

        shackle_line = Rectangle(
            width=shackle_width, height=shackle_height-shackle_radius,
            stroke_width=0, fill_opacity=1
        ).next_to(body, UP, buff=0)

        shackle_arc = AnnularSector(
            inner_radius=shackle_radius - shackle_width/2,
            outer_radius=shackle_radius + shackle_width/2,
            angle=PI,
        ).next_to(shackle_line, UP, buff=0)

        super().__init__(
            body, shackle_arc,
            shackle_line.copy().shift(shackle_radius*LEFT),
            shackle_line.copy().shift(shackle_radius*RIGHT),
            **kwargs
        )

        self.set(color=color)


class Pointers(VMobject):
    def __init__(self, num_pointers=2, width=1.2, height=1, buff=0, **kwargs):
        super().__init__(**kwargs)
        pointer_height = (height - (num_pointers - 1) * buff) / num_pointers
        for i in range(num_pointers):
            base = i * (pointer_height + buff) * UP
            self.start_new_path(base)
            self.add_line_to(base + width/2 * RIGHT + pointer_height * UP)
            self.add_line_to(base + width * RIGHT)


class IdealGasLawScene(Scene):
    def construct(self):
        N = 256
        T = ValueTracker(1)

        time_scale = ValueTracker(1)

        container = Rectangle(width=8, height=6, stroke_color=WHITE, fill_opacity=0.05, fill_color=NORMAL_TEMP_COLOR)

        def get_color():
            if T.get_value() < 1:
                return interpolate_color(LOW_TEMP_COLOR, NORMAL_TEMP_COLOR, 2*sigmoid(np.log(T.get_value())))
            else:
                return interpolate_color(NORMAL_TEMP_COLOR, HIGH_TEMP_COLOR, 2*sigmoid(np.log(T.get_value())) - 1)
        
        container.add_updater(lambda mob: mob.set(fill_color=get_color()), call_updater=True)

        self.add(container)

        particle = Dot(radius=0.04, color=NORMAL_TEMP_COLOR)

        size_diff = container.stroke_width/200 + particle.radius

        show_pressure = False

        def add_pressure_vec(x, y, dir, size):
            pressure_vec = Vector(0.5*size*dir, color=PRESSURE_COLOR)
            pressure_vec.shift(x*RIGHT + y*UP + (size_diff + container.stroke_width/200)*dir)

            t = 0
            def pressure_vec_updater(mob, dt):
                nonlocal t
                t += time_scale.get_value() * dt
                if t >= 1:
                    self.remove(mob)
                else:
                    mob.set_opacity((1 - t)**3)

            pressure_vec.add_updater(pressure_vec_updater)

            self.add(pressure_vec)

        def particle_updater():
            p = np.random.uniform(-1, 1) * container.width/2 * RIGHT + np.random.uniform(-1, 1) * container.height/2 * UP 
            base_v = np.random.uniform(-1, 1) * RIGHT + np.random.uniform(-1, 1) * UP

            def updater(mob, dt):
                nonlocal p, base_v

                v = base_v * np.sqrt(T.get_value())
                p += v * time_scale.get_value() * dt 

                if p[0] < (x := container.get_left()[0] + size_diff):
                    if show_pressure:
                        y = p[1] - v[1] * (p[0] - x) / v[0]
                        add_pressure_vec(x, y, LEFT, np.sqrt(abs(v[0])))
                    p[0] = 2*x - p[0]
                    base_v[0] *= -1
                if p[0] > (x := container.get_right()[0] - size_diff):
                    if show_pressure:
                        y = p[1] - v[1] * (p[0] - x) / v[0]
                        add_pressure_vec(x, y, RIGHT, np.sqrt(abs(v[0])))
                    p[0] = 2*x - p[0]
                    base_v[0] *= -1
                if p[1] < (y := container.get_bottom()[1] + size_diff):
                    if show_pressure:
                        x = p[0] - v[0] * (p[1] - y) / v[1]
                        add_pressure_vec(x, y, DOWN, np.sqrt(abs(v[1])))
                    p[1] = 2*y - p[1]
                    base_v[1] *= -1
                if p[1] > (y := container.get_top()[1] - size_diff):
                    if show_pressure:
                        x = p[0] - v[0] * (p[1] - y) / v[1]
                        add_pressure_vec(x, y, UP, np.sqrt(abs(v[1])))
                    p[1] = 2*y - p[1]
                    base_v[1] *= -1
                    
                mob.move_to(p)

                mob.set(color=get_color())

            return updater

        particles = Group(*[particle.copy().add_updater(particle_updater()) for _ in range(N)])
            
        self.add(particles)

        self.wait(5)

        container_buff = (self.camera.frame_height - container.height) / 2

        self.play(
            container.animate.to_edge(LEFT, buff=container_buff),
            run_time=4
        )

        self.wait(2)

        box_width = self.camera.frame_width - container.width - 3*container_buff
        box_buff = 0.12

        box_config = {
            'height': 0.8,
            'corner_radius': 0.025,
            'fill_opacity': 0.05,
            'stroke_width': 1,
        }

        sub_box_width = (box_width - box_buff) / 2

        gas_eq_box = RoundedRectangle(width=box_width, **box_config)
        P_box = RoundedRectangle(width=sub_box_width, color=PRESSURE_COLOR, **box_config)
        V_box = RoundedRectangle(width=sub_box_width, color=VOLUME_COLOR, **box_config)
        N_box = RoundedRectangle(width=sub_box_width, color=NORMAL_TEMP_COLOR, **box_config)
        T_box = RoundedRectangle(width=sub_box_width, color=NORMAL_TEMP_COLOR, **box_config)

        T_box.add_updater(lambda mob: mob.set(color=get_color()))

        Group(
            gas_eq_box,
            Group(
                P_box, V_box, N_box, T_box
            ).arrange_in_grid(2, 2, buff=box_buff)
        ).arrange(DOWN, buff=2*box_buff).to_edge(RIGHT, buff=container_buff)

        gas_eq = MathTex('P', 'V', '=', 'N', 'k_B', 'T', background_stroke_width=4).scale(0.8)
        gas_eq[0].set(color=PRESSURE_COLOR)
        gas_eq[1].set(color=VOLUME_COLOR)
        gas_eq[3].set(color=NORMAL_TEMP_COLOR)
        gas_eq[5].add_updater(lambda mob: mob.set(color=get_color()), call_updater=True)

        gas_eq.move_to(gas_eq_box)

        P_label = gas_eq[0].copy()
        V_label = gas_eq[1].copy()
        N_label = gas_eq[3].copy()
        T_label = gas_eq[5].copy()

        self.play(
            FadeIn(Group(gas_eq, gas_eq_box), shift=0.5*DOWN),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(3)

        self.play(AnimationGroup(
            FadeIn(P_box, shift=box_buff*DOWN),
            P_label.animate.move_to(P_box),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        show_pressure = True

        self.wait(3)

        self.play(AnimationGroup(
            FadeIn(V_box, shift=box_buff*DOWN),
            V_label.animate.move_to(V_box),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        container_trace = container.copy().set(color=VOLUME_COLOR).set_fill(opacity=0)

        self.play(
            Create(container_trace),
            run_time=3, rate_func=rate_functions.ease_in_out_cubic
        )

        self.remove(container_trace)
        container.set_stroke(color=VOLUME_COLOR)

        self.wait(3)

        self.play(AnimationGroup(
            FadeIn(N_box, shift=box_buff*DOWN, rate_func=rate_functions.ease_out_cubic),
            N_label.animate(rate_func=rate_functions.ease_out_cubic).move_to(N_box),
            time_scale.animate(rate_func=rate_functions.ease_in_quad).set_value(0),
            lag_ratio=0.1, run_time=3
        ))

        self.wait()

        circles = [Circle(radius=0.1, stroke_width=1, color=NORMAL_TEMP_COLOR).move_to(p) for p in particles]

        self.play(AnimationGroup(
            *[Create(circle) for circle in circles],
            lag_ratio=2, run_time=5, rate_func=rate_functions.ease_in_quad
        ))
        
        self.wait(2)

        self.play(AnimationGroup(
            *[FadeOut(circle) for circle in circles],
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.play(
            time_scale.animate.set_value(1),
            run_time=4, rate_func=rate_functions.ease_in_quad
        )

        self.wait()

        self.play(AnimationGroup(
            FadeIn(T_box, shift=box_buff*DOWN),
            T_label.animate.move_to(T_box),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        lock = Lock().scale_to_fit_height(P_label.height)

        P_lock = lock.copy().set(color=PRESSURE_COLOR).move_to(P_box)
        V_lock = lock.copy().set(color=VOLUME_COLOR).move_to(V_box)
        N_lock = lock.copy().set(color=NORMAL_TEMP_COLOR).move_to(N_box)
        T_lock = lock.copy().set(color=NORMAL_TEMP_COLOR).move_to(T_box)

        shift_width = 0.2 * sub_box_width

        Group(P_lock, V_lock, N_lock, T_lock).shift(shift_width*RIGHT)
        T_lock.add_updater(lambda mob: mob.set(color=get_color()))

        P_colon = MathTex(':', color=PRESSURE_COLOR).scale(0.8).move_to(P_box)
        V_colon = MathTex(':', color=VOLUME_COLOR).scale(0.8).move_to(V_box)
        N_colon = MathTex(':', color=NORMAL_TEMP_COLOR).scale(0.8).move_to(N_box)
        T_colon = MathTex(':', color=NORMAL_TEMP_COLOR).scale(0.8).move_to(T_box)

        T_colon.add_updater(lambda mob: mob.set(color=get_color()))

        self.play(AnimationGroup(
            AnimationGroup(
                N_label.animate.shift(shift_width*LEFT),
                FadeIn(N_colon, shift=shift_width*LEFT)
            ),
            FadeIn(N_lock, scale=1.5),
            AnimationGroup(
                V_label.animate.shift(shift_width*LEFT),
                FadeIn(V_colon, shift=shift_width*LEFT)
            ),
            AnimationGroup(
                FadeIn(V_lock, scale=1.5),
                container.animate.set_stroke(color=WHITE)
            ),
            lag_ratio=0.1, run_time=5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        up = Pointers(buff=-0.1, stroke_width=10*lock.height).scale_to_fit_height(lock.height)

        P_up = up.copy().set(color=PRESSURE_COLOR).move_to(P_lock)
        V_up = up.copy().set(color=VOLUME_COLOR).move_to(V_lock)
        N_up = up.copy().set(color=NORMAL_TEMP_COLOR).move_to(N_lock)
        T_up = up.copy().set(color=NORMAL_TEMP_COLOR).move_to(T_lock)

        T_up.add_updater(lambda mob: mob.set(color=get_color()))

        self.play(AnimationGroup(
            AnimationGroup(
                AnimationGroup(
                    T_label.animate.shift(shift_width*LEFT),
                    FadeIn(T_colon, shift=shift_width*LEFT)
                ),
                FadeIn(T_up, shift=up.height*UP),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            T.animate(rate_func=rate_functions.ease_in_cubic).set_value(5),
            lag_ratio=0.3, run_time=6
        ))

        v_eq = MathTex(r'(v \propto \sqrt{T})', color=GRAY).scale(0.67).next_to(T_box, DOWN, buff=SMALL_BUFF)

        self.play(
            FadeIn(v_eq, target_position=T_box.get_bottom()),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        down = up.copy().rotate(180*DEGREES)

        P_down = down.copy().set(color=PRESSURE_COLOR).move_to(P_lock)
        V_down = down.copy().set(color=VOLUME_COLOR).move_to(V_lock)
        N_down = down.copy().set(color=NORMAL_TEMP_COLOR).move_to(N_lock)
        T_down = down.copy().set(color=NORMAL_TEMP_COLOR).move_to(T_lock)

        T_down.add_updater(lambda mob: mob.set(color=get_color()))

        self.play(AnimationGroup(
            AnimationGroup(
                P_label.animate.shift(shift_width*LEFT),
                FadeIn(P_colon, shift=shift_width*LEFT)
            ),
            FadeIn(P_up, shift=up.height*UP),
            lag_ratio=0.2, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(5)

        dot = Dot(radius=0.05)

        P_dot = dot.copy().set(color=PRESSURE_COLOR).move_to(P_lock)
        V_dot = dot.copy().set(color=VOLUME_COLOR).move_to(V_lock)
        N_dot = dot.copy().set(color=NORMAL_TEMP_COLOR).move_to(N_lock)
        T_dot = dot.copy().set(color=NORMAL_TEMP_COLOR).move_to(T_lock)

        T_dot.add_updater(lambda mob: mob.set(color=get_color()))

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(T_up, shift=up.height*DOWN),
                FadeIn(T_dot, shift=up.height*DOWN),
                rate_func=rate_functions.ease_out_cubic
            ),
            T.animate(rate_func=rate_functions.linear).set_value(1),
            AnimationGroup(
                FadeOut(P_up, shift=up.height*DOWN),
                FadeIn(P_dot, shift=up.height*DOWN),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=5
        ))

        self.play(AnimationGroup(
            AnimationGroup(
                FadeIn(T_down, shift=down.height*DOWN),
                FadeOut(T_dot, shift=down.height*DOWN),
            ),
            T.animate(rate_func=rate_functions.ease_in_cubic).set_value(0.2),
            lag_ratio=0.3, run_time=6,
        ))

        self.play(AnimationGroup(
            FadeIn(P_down, shift=down.height*DOWN),
            FadeOut(P_dot, shift=down.height*DOWN),
            run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(5)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(T_down, shift=down.height*UP),
                FadeIn(T_dot, shift=down.height*UP),
                rate_func=rate_functions.ease_out_cubic
            ),
            T.animate(rate_func=rate_functions.linear).set_value(1),
            AnimationGroup(
                FadeOut(P_down, shift=down.height*UP),
                FadeIn(P_dot, shift=down.height*UP),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=5
        ))

        self.wait(3)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(T_dot, scale=0.67),
                FadeIn(T_lock, scale=1.5),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                AnimationGroup(
                    FadeOut(V_lock, scale=1.5),
                    container.animate.set_stroke(color=VOLUME_COLOR),
                ),
                FadeIn(V_dot, scale=0.67),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.05, run_time=5
        ))

        self.wait(3)

        show_pressure = False

        self.play(AnimationGroup(
            AnimationGroup(
                FadeIn(V_down, shift=down.height*DOWN),
                FadeOut(V_dot, shift=down.height*DOWN),
                rate_func=rate_functions.ease_out_cubic
            ),
            container.animate(rate_func=rate_functions.ease_in_out_cubic).scale(0.5),
            lag_ratio=0.2, run_time=6,
        ))

        show_pressure = True

        self.play(AnimationGroup(
            FadeIn(P_up, shift=up.height*UP),
            FadeOut(P_dot, shift=up.height*UP),
            run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(5)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(V_down, shift=down.height*UP),
                FadeIn(V_dot, shift=down.height*UP),
                rate_func=rate_functions.ease_out_cubic
            ),
            container.animate(rate_func=rate_functions.ease_in_out_cubic).scale(2),
            AnimationGroup(
                FadeOut(P_up, shift=up.height*DOWN),
                FadeIn(P_dot, shift=up.height*DOWN),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=5
        ))

        self.wait(3)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(P_dot, scale=0.67),
                FadeIn(P_lock, scale=1.5),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeOut(N_lock, scale=1.5),
                FadeIn(N_dot, scale=0.67),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.05, run_time=5
        ))

        self.wait(3)

        show_pressure = False

        self.play(AnimationGroup(
            AnimationGroup(
                *[FadeOut(p) for p in particles[::2]],
                lag_ratio=0.005, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(N_down, shift=down.height*DOWN),
                FadeOut(N_dot, shift=down.height*DOWN),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(V_down, shift=up.height*UP),
                FadeOut(V_dot, shift=up.height*UP),
                rate_func=rate_functions.ease_out_cubic
            ),
            container.animate(rate_func=rate_functions.ease_out_cubic).scale(1/np.sqrt(2)),
            lag_ratio=0.2, run_time=7.5,
        ))

        show_pressure = True

        self.wait(6)
