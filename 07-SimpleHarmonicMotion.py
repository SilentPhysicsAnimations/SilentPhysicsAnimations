from manim import *

#########################

DISPLACEMENT_COLOR = BLUE

#########################

class Spring(ParametricFunction):
    def __init__(self, length=1, radius=0.1, coils=4, **kwargs):
        self.length = length
        self.radius = radius

        self._center_line_length = length - 2*radius

        t_max = 2 * PI * (coils + 0.5)

        super().__init__(
            function=lambda t: np.array([
                radius * np.cos(t+PI) + self._center_line_length * (t / t_max),
                radius * np.sin(t+PI),
                0
            ]),
            t_range=[0, t_max],
            **kwargs
        )

    def compress(self, x):
        self.stretch_to_fit_width(self.length - x, about_point=self.get_start())
        return self

    def get_compression(self):
        return self.length - self.width
        

class SimpleHarmonicMotionScene(MovingCameraScene):
    def construct(self):
        ball = Dot(radius=0.5, z_index=1)
        spring = Spring(length=4, radius=0.2, coils=8, z_index=1).next_to(ball, LEFT, buff=0)
        ball.add_updater(lambda mob: mob.next_to(spring, buff=0))

        wall = Line(ORIGIN, DOWN).next_to(spring, LEFT, buff=0)

        self.add(spring, ball, wall)

        A = 2
        k = 20
        m = 3
        omega = np.sqrt(k/m)
        T = 2*PI / omega

        def x(t):
            return A * np.cos(omega*t)
        
        time = ValueTracker(-3*T/4)
        def spring_updater(mob):
            t = time.get_value()
            mob.compress(-x(t))
            
        spring.add_updater(spring_updater)

        center_line = DashedLine(10*UP, 10*DOWN, stroke_width=1, color=GRAY)
        left_line = center_line.copy().shift(A*LEFT)
        right_line = center_line.copy().shift(A*RIGHT)

        self.wait()

        self.play(
            Create(center_line),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )
        
        self.wait()

        spring_brace = BraceBetweenPoints(
            ball.get_bottom() + A*LEFT,
            ball.get_bottom(),
            color=DISPLACEMENT_COLOR
        )
        spring_brace.stretch_to_fit_width(1e-6, about_point=spring_brace.get_right())

        self.play(AnimationGroup(
            time.animate(rate_func=rate_functions.linear)
                .increment_value(T/4),
            spring_brace.animate(rate_func=rate_functions.ease_out_sine)
                .stretch_to_fit_width(A, about_point=spring_brace.get_right()),
            run_time=1
        ))

        self.wait(0.2)

        spring_brace_label = MathTex('A', color=DISPLACEMENT_COLOR).next_to(spring_brace, DOWN)

        self.play(
            FadeIn(spring_brace_label, target_position=spring_brace),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(
            Create(left_line),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(
            time.animate.increment_value(T),
            run_time=T,
            rate_func=rate_functions.linear
        )

        self.play(
            time.animate.increment_value(T/2),
            run_time=(T/2) * (PI/2),
            rate_func=rate_functions.ease_out_sine
        )

        self.wait(0.5)

        spring_brace_copy = spring_brace.copy()
        spring_brace_label_copy = spring_brace_label.copy()

        self.add(spring_brace_copy, spring_brace_label_copy)

        self.play(AnimationGroup(
            AnimationGroup(
                Group(spring_brace_copy, spring_brace_label_copy).animate.shift(A*RIGHT),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                Create(right_line),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=3.5
        ))

        self.wait()

        time.set_value(0)

        trace = VMobject()
        trace.add_updater(lambda mob: mob.become(
            ParametricFunction(
                lambda t: [x(t), 0.8*(time.get_value() - t) + 0.5, 0], t_range=[0, time.get_value()],
                stroke_width=2, color=DISPLACEMENT_COLOR
            )
        ))

        self.add(trace)
        time.add_updater(lambda vt, dt: vt.increment_value(dt))

        self.play(
            self.camera.frame.animate.scale(1.5).shift(2.5*UP),
            run_time=3,
        )

        self.wait_until(lambda: time.get_value() >= 5*T)
        time.clear_updaters()
        time.set_value(5*T)

        self.play(
            time.animate.increment_value(T/2),
            run_time=(T/2) * (PI/2),
            rate_func=rate_functions.ease_out_sine
        )

        trace.clear_updaters()

        self.wait(0.5)

        axes = Axes(
            x_range=[0, 5.5*T],
            y_range=[-1.5*A, 1.5*A],
            x_length=5.5,
            y_length=4,
            axis_config={
                'include_ticks': False,
                'tip_width': 0.2,
                'tip_height': 0.2,
            }
        )

        axes.next_to(right_line, buff=3)

        t_label = MathTex('t').next_to(axes.get_x_axis(), RIGHT)
        x_label = MathTex('x').next_to(axes.get_y_axis(), UP)

        axes.add(t_label, x_label)

        real_curve = axes.plot(x, stroke_width=2, color=DISPLACEMENT_COLOR)

        self.add(axes)

        self.play(AnimationGroup(
            self.camera.frame.animate.move_to(Group(wall, axes)),
            FadeIn(axes, shift=3*LEFT),
            run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(
            ReplacementTransform(trace.copy(), real_curve),
            run_time=3
        )
        
        self.wait(2)

        self.play(AnimationGroup(
            self.camera.frame.animate.move_to(axes).scale(0.75),
            FadeOut(Group(
                spring, wall, ball, right_line, center_line, left_line, trace,
                spring_brace, spring_brace_copy, spring_brace_label, spring_brace_label_copy
            )),
            run_time=2.5
        ))

        self.wait(2)

        guess_curve = axes.plot(np.cos)
        guess_cos = MathTex(r'\cos', '(', 't', ')').next_to(axes, MED_LARGE_BUFF*DOWN)

        self.play(AnimationGroup(
            Create(guess_curve, rate_func=rate_functions.ease_out_cubic),
            FadeIn(guess_cos, shift=DOWN, rate_func=rate_functions.ease_out_cubic),
            lag_ratio=0.5, run_time=3
        ))

        self.wait(1.5)

        guess_A = MathTex('A', r'\cos', '(', 't', ')').move_to(guess_cos)

        self.play(AnimationGroup(
            guess_curve.animate.stretch_to_fit_height(real_curve.height),
            AnimationGroup(
                FadeIn(guess_A[0], shift=0.1*LEFT),
                guess_cos.animate.move_to(guess_A[1:])
            ),
            lag_ratio=0.2, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*guess_cos)
        self.add(guess_A)

        self.wait()

        guess_omega = MathTex('A', r'\cos', '(', r'\omega', 't', ')').move_to(guess_cos)

        omega_tracker = ValueTracker(1)
        guess_curve.add_updater(lambda mob: mob.become(
            axes.plot(lambda t: A*np.cos(omega_tracker.get_value()*t))
        ))

        self.play(AnimationGroup(
            omega_tracker.animate.set_value(omega),
            AnimationGroup(
                FadeIn(guess_omega[3], shift=0.1*DOWN),
                guess_A[0:3].animate.move_to(guess_omega[0:3]),
                guess_A[3:].animate.move_to(guess_omega[4:])
            ),
            lag_ratio=0.2, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*guess_A)
        self.add(guess_omega)

        guess_curve.clear_updaters()

        self.wait()

        final_eq = MathTex('x', '(', 't', ')', '=', 'A', r'\cos', '(', r'\omega', 't', ')', color=DISPLACEMENT_COLOR).move_to(guess_cos)

        self.play(AnimationGroup(
            FadeOut(guess_curve),
            AnimationGroup(
                FadeIn(final_eq[:5], shift=0.3*LEFT),
                guess_omega.animate.move_to(final_eq[5:]).set(color=DISPLACEMENT_COLOR)
            ),
            lag_ratio=0.5, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*guess_omega)
        self.add(final_eq)

        self.wait(3)

        time.add_updater(lambda tracker, dt: tracker.increment_value(dt))

        self.play(AnimationGroup(
            self.camera.frame.animate.move_to(Group(wall, axes)).scale(1.2),
            FadeIn(Group(spring, wall, ball)),
            run_time=2.5
        ))

        self.wait(2)

        k_line = NumberLine(x_range=[0, 2*k], length=3, include_ticks=False).shift(4*DOWN)
        k_line.add(k_line.get_tick(0), k_line.get_tick(k, size=0.05), k_line.get_tick(2*k))
        k_line_dot = Dot(radius=0.075).move_to(k_line)
        k_line_label = MathTex('k').next_to(k_line, LEFT)

        m_line = NumberLine(x_range=[0, 2*m], length=3, include_ticks=False).shift(5*DOWN)
        m_line.add(m_line.get_tick(0), m_line.get_tick(m, size=0.05), m_line.get_tick(2*m))
        m_line_dot = Dot(radius=0.075).move_to(m_line)
        m_line_label = MathTex('m').next_to(m_line, LEFT)

        omega_line = NumberLine(x_range=[0, 2*omega], length=3, include_ticks=False).shift(4.5*DOWN).set_x(axes.get_x())
        omega_line.add(omega_line.get_tick(0), omega_line.get_tick(omega, size=0.05), omega_line.get_tick(2*omega))
        omega_line_dot = Dot(radius=0.075).move_to(omega_line)
        omega_line_label = MathTex(r'\omega').next_to(omega_line, LEFT)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(1.5*DOWN),
            FadeIn(Group(k_line, k_line_dot, k_line_label), shift=0.5*DOWN),
            FadeIn(Group(m_line, m_line_dot, m_line_label), shift=0.5*DOWN),
            lag_ratio=0.1, run_time=2.5
        ))

        self.wait(0.5)

        self.play(
            FadeIn(Group(omega_line, omega_line_dot, omega_line_label), shift=0.5*RIGHT),
            run_time=1.5
        )

        k_tracker = ValueTracker(k)
        m_tracker = ValueTracker(m)

        k_line_dot.add_updater(lambda mob: mob.move_to(k_line.n2p(k_tracker.get_value())))
        m_line_dot.add_updater(lambda mob: mob.move_to(m_line.n2p(m_tracker.get_value())))
        omega_line_dot.add_updater(lambda mob: mob.move_to(omega_line.n2p(
            np.sqrt(k_tracker.get_value() / m_tracker.get_value())
        )))

        spring.add_updater(lambda mob: mob.set(stroke_width=4 * (k_tracker.get_value() / k)))
        ball.add_updater(lambda mob: mob.scale_to_fit_width(m_tracker.get_value() / m))

        time.clear_updaters()
        time.add_updater(lambda tracker, dt: tracker.increment_value(
            (np.sqrt(k_tracker.get_value() / m_tracker.get_value()) / omega) * dt
        ))

        real_curve.add_updater(lambda mob: mob.become(axes.plot(
            lambda t: A * np.cos(np.sqrt(k_tracker.get_value() / m_tracker.get_value()) * t),
            stroke_width=2, color=DISPLACEMENT_COLOR
        )))

        self.wait()

        k_line_dot.save_state()
        self.play(
            k_line_dot.animate.scale(1.5),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.2)

        self.play(
            k_tracker.animate.set_value(1.8*k),
            omega_tracker.animate.set_value(np.sqrt(1.8)*omega),
            run_time=3
        )
        self.play(
            k_tracker.animate.set_value(0.2*k),
            omega_tracker.animate.set_value(np.sqrt(0.2)*omega),
            run_time=3
        )
        self.play(
            k_tracker.animate.set_value(k),
            omega_tracker.animate.set_value(omega),
            run_time=3
        )

        self.play(
            Restore(k_line_dot),
            rate_func=rate_functions.ease_out_cubic
        )

        omega_prop_sqrt_k = MathTex(r'\omega', r'\propto', r'\sqrt{k}')
        omega_prop_sqrt_k.shift(k_line.get_right() + 0.5*RIGHT - omega_prop_sqrt_k[0].get_left())

        self.play(
            FadeIn(omega_prop_sqrt_k, shift=0.25*RIGHT),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        m_line_dot.save_state()
        self.play(
            m_line_dot.animate.scale(1.5),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.2)

        self.play(
            m_tracker.animate.set_value(1.8*m),
            omega_tracker.animate.set_value(np.sqrt(1.8)*omega),
            run_time=3
        )
        self.play(
            m_tracker.animate.set_value(0.2*m),
            omega_tracker.animate.set_value(np.sqrt(0.2)*omega),
            run_time=3
        )
        self.play(
            m_tracker.animate.set_value(m),
            omega_tracker.animate.set_value(omega),
            run_time=3
        )

        self.play(
            Restore(m_line_dot),
            rate_func=rate_functions.ease_out_cubic
        )

        omega_prop_inv_sqrt_m = MathTex(r'\omega', r'\propto', r'\frac{1}{\sqrt{m}}')
        omega_prop_inv_sqrt_m.shift(m_line.get_right() + 0.5*RIGHT - omega_prop_inv_sqrt_m[0].get_left())

        self.play(
            FadeIn(omega_prop_inv_sqrt_m, shift=0.25*RIGHT),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        omega_eq = MathTex(r'\omega', '=', r'\sqrt{\frac{k}{m}}').move_to(omega_line)

        self.play(AnimationGroup(
            FadeOut(Group(omega_line, omega_line_dot), shift=0.5*RIGHT),
            omega_line_label.animate.move_to(omega_eq[0]),
            FadeIn(omega_eq[1:], shift=0.2*RIGHT),
            lag_ratio=0.2, run_time=3.5
        ))

        self.remove(omega_line_label)
        self.add(omega_eq)

        while time.get_value() > 0:
            time.increment_value(-T)

        self.wait_until(lambda: time.get_value() >= 0)
        time.clear_updaters()
        time.set_value(0)

        self.play(
            time.animate.increment_value(5*T/4),
            run_time=(5*T/4) * (PI/2),
            rate_func=rate_functions.ease_out_sine
        )

        self.wait()
