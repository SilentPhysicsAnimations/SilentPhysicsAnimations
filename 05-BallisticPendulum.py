from manim import *

########################

DISPLACEMENT_COLOR = BLUE
VELOCITY_COLOR = RED
EQUALITY_COLOR = YELLOW

########################

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


class PendulumHead(Arc):
    def __init__(
        self,
        inner_radius=1,
        outer_radius=2,
        angle=TAU / 4,
        start_angle=0,
        fill_opacity=1,
        stroke_width=0,
        color=WHITE,
        **kwargs,
    ):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        super().__init__(
            start_angle=start_angle,
            angle=angle,
            fill_opacity=fill_opacity,
            stroke_width=stroke_width,
            color=color,
            **kwargs,
        )

    def generate_points(self):
        inner_arc = Arc(
            start_angle=-PI/2,
            angle=PI,
            radius=self.inner_radius,
            arc_center=self.arc_center
        )
        angle_delta = np.arccos(self.inner_radius/self.outer_radius)
        outer_arc = Arc(
            start_angle=-PI/2 - angle_delta,
            angle=PI + 2*angle_delta,
            radius=self.outer_radius,
            arc_center=self.arc_center
        )
        outer_arc.reverse_points()
        self.append_points(inner_arc.points)
        self.add_line_to(outer_arc.points[0])
        self.append_points(outer_arc.points)
        self.add_line_to(inner_arc.points[0])

    init_points = generate_points


class BallisticPendulumScene(MovingCameraScene):
    def construct(self):
        fill = {
            'fill_color': WHITE,
            'fill_opacity': 1,
            'stroke_width': 0,
        }

        k = 300
        spring = Spring(length=2, coils=8).shift(3*LEFT + 2*DOWN)

        inner_radius = 0.25
        outer_radius = 0.45

        m = 1
        bullet = Circle(radius=inner_radius + 1e-3, z_index=0, **fill).next_to(spring, buff=0)
        bullet.set_sheen(-0.1, LEFT)
        
        l = 2.5
        M = 3
        pendulum_back = Circle(radius=outer_radius, stroke_width=0, fill_color=WHITE, fill_opacity=0.9)
        pendulum_back.set_sheen(-0.5, RIGHT)
        pendulum_back.set_y(bullet.get_y())
        pendulum_head = PendulumHead(inner_radius=inner_radius, outer_radius=outer_radius)
        pendulum_head.set_sheen(-0.1, RIGHT)
        pendulum_head.shift(pendulum_back.get_center() - pendulum_head.arc_center),
        pendulum_shaft = Line(ORIGIN, l*DOWN).next_to(pendulum_back, UP, buff=0)
        pendulum = Group(pendulum_back, pendulum_head, pendulum_shaft).shift(3*RIGHT)

        self.camera.frame.move_to(Group(spring, bullet, pendulum)).scale(0.9)

        self.add(spring, bullet, pendulum)

        self.wait()

        delta_x = 1

        bullet.add_updater(lambda mob: mob.next_to(spring, buff=0))

        self.play(
            spring.animate.compress(delta_x),
            run_time=3, rate_func=rate_functions.ease_out_quad
        )

        self.wait()

        omega_spring = np.sqrt(k/m)

        self.play(
            spring.animate.compress(0),
            run_time=PI/(2*omega_spring), rate_func=rate_functions.ease_in_sine
        )

        dist = pendulum_back.get_x() - bullet.get_x()
        v = delta_x * omega_spring
        bullet.clear_updaters()

        self.play(
            bullet.animate.shift(dist*RIGHT),
            run_time=dist/v, rate_func=rate_functions.linear
        )

        pendulum.add(bullet)

        g = 9.82
        u = m*v / (m + M)
        h = u*u / (2*g)

        delta_theta = np.arccos(1 - h/l)
        omega_pendulum = np.sqrt(g/l)

        self.play(
            Rotate(pendulum, delta_theta, about_point=pendulum_shaft.get_start()),
            run_time=PI/(2*omega_pendulum), rate_func=rate_functions.ease_out_sine
        )

        self.play(
            Rotate(pendulum, -delta_theta, about_point=pendulum_shaft.get_start()),
            run_time=PI/(2*omega_pendulum), rate_func=rate_functions.ease_in_sine
        )

        self.wait(2)

        pendulum.remove(bullet)

        self.play(AnimationGroup(
            bullet.animate.next_to(spring, buff=0),
            self.camera.frame.animate.move_to(spring.get_right()).shift(0.5*delta_x*LEFT).scale(5/9),
            lag_ratio=0.1, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        spring_corner = bullet.get_bottom()
        spring_brace = BraceBetweenPoints(
            spring_corner + 2*delta_x*LEFT,
            spring_corner,
            buff=SMALL_BUFF,
            color=DISPLACEMENT_COLOR
        )
        spring_brace.scale(0.5, about_point=spring_brace.get_corner(UR))
        spring_brace.stretch_to_fit_width(1e-6, about_point=spring_brace.get_right())

        self.play(
            FadeIn(spring_brace, target_position=bullet),
            run_time=0.5, rate_func=rate_functions.ease_out_cubic
        )

        bullet.add_updater(lambda mob: mob.next_to(spring, buff=0))

        self.play(
            spring.animate.compress(delta_x),
            spring_brace.animate.stretch_to_fit_width(delta_x, about_point=spring_brace.get_right()),
            run_time=3, rate_func=rate_functions.ease_out_quad
        )

        self.wait(0.3)

        spring_brace_label = MathTex('x', color=DISPLACEMENT_COLOR).scale(0.4).next_to(spring_brace, DOWN, buff=SMALL_BUFF)
        spring_constant_label = MathTex('k').scale(0.4).next_to(spring.get_corner(UL), UP, buff=SMALL_BUFF)
        
        self.play(AnimationGroup(
            FadeIn(spring_brace_label, target_position=spring_brace),
            FadeIn(spring_constant_label, target_position=spring.get_left()),
            lag_ratio=0.5, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        circled_one = Group(
            Circle(radius=inner_radius, stroke_width=1, color=WHITE),
            MathTex('1')
        ).scale(0.5).move_to(bullet)

        bullet.set(z_index=1)
        self.play(
            circled_one.animate.shift(MED_LARGE_BUFF*UP),
            rate_func=rate_functions.ease_out_cubic
        )
        bullet.set(z_index=0)

        self.wait(2)

        energy_eq_1 = MathTex('E_1', '=', r'\tfrac{1}{2}', 'k', 'x', '^2')
        energy_eq_1[4].set(color=DISPLACEMENT_COLOR)
        energy_eq_1.scale(0.4).next_to(circled_one, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            FadeIn(energy_eq_1[:3], shift=SMALL_BUFF*RIGHT),
            AnimationGroup(
                spring_constant_label.animate.move_to(energy_eq_1[3]),
                spring_brace_label.animate.move_to(energy_eq_1[4]),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            FadeIn(energy_eq_1[5], target_position=energy_eq_1[4]),
            lag_ratio=0.2, run_time=3.5
        ))

        self.add(energy_eq_1)
        self.remove(spring_constant_label, spring_brace_label)
        
        self.wait(2)

        vector_config = {
            'max_tip_length_to_length_ratio': 0.15,
            'stroke_width': 3,
            'max_stroke_width_to_length_ratio': 1000
        }

        v_2_vector = Vector(color=VELOCITY_COLOR, **vector_config).shift(999*RIGHT)

        vector_scale = 0.08

        def v_2_vector_updater(mob):
            x = spring.get_compression()
            v = omega_spring * np.sqrt(delta_x**2 - x**2) + 1e-6
            mob.put_start_and_end_on(ORIGIN, vector_scale * v * RIGHT)
            mob.shift(bullet.get_right())

        v_2_vector.add_updater(v_2_vector_updater)

        self.add(v_2_vector)

        slow_down_factor = 8

        self.play(
            spring.animate.compress(0),
            run_time=slow_down_factor*PI/(2*omega_spring), rate_func=rate_functions.ease_in_sine
        )

        bullet.clear_updaters()

        bullet_speed = ValueTracker(v/slow_down_factor)

        def slow_down_updater(mob, dt):
            mob.shift(bullet_speed.get_value() * dt * RIGHT)

        bullet.add_updater(slow_down_updater)

        self.play(
            self.camera.frame.animate.scale(0.9).move_to(bullet),
            bullet_speed.animate.set_value(0.067),
            run_time=0.8, rate_func=rate_functions.ease_out_cubic
        )

        v_2_vector_label = MathTex('v', '_2', color=VELOCITY_COLOR).scale(0.4).next_to(v_2_vector, RIGHT, buff=SMALL_BUFF)
        v_2_vector_label.add_updater(slow_down_updater)
        mass_label = MathTex('m', color=BLACK).scale(0.4).add_updater(lambda mob: mob.move_to(bullet))

        self.play(AnimationGroup(
            FadeIn(v_2_vector_label),
            FadeIn(mass_label),
            lag_ratio=0.5, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        circled_two = Group(
            Circle(radius=inner_radius, stroke_width=1, color=WHITE),
            MathTex('2')
        ).scale(0.5).move_to(bullet)

        self.play(
            circled_two.animate.shift(MED_LARGE_BUFF*UP),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        energy_eq_2 = MathTex('E_2', '=', r'\tfrac{1}{2}', 'm', 'v', '^2', '_2')
        energy_eq_2[4::2].set(color=VELOCITY_COLOR)
        energy_eq_2.scale(0.4).next_to(circled_two, buff=SMALL_BUFF)

        mass_label.clear_updaters()
        v_2_vector_label.clear_updaters()

        self.play(AnimationGroup(
            FadeIn(energy_eq_2[:3], shift=SMALL_BUFF*RIGHT),
            AnimationGroup(
                mass_label.animate.move_to(energy_eq_2[3]).set(color=WHITE),
                AnimationGroup(
                    v_2_vector_label[0].animate.move_to(energy_eq_2[4]),
                    v_2_vector_label[1].animate.move_to(energy_eq_2[6]),
                ),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            FadeIn(energy_eq_2[5], target_position=energy_eq_2[4::2]),
            lag_ratio=0.15, run_time=3.5
        ))

        self.add(energy_eq_2)
        self.remove(mass_label, *v_2_vector_label)
        
        self.wait(2)

        pendulum.shift(4*RIGHT)

        slow_down_factor = 4

        self.play(
            self.camera.frame.animate(rate_func=rate_functions.ease_in_out_cubic, run_time=2).scale(1.2).move_to(pendulum),
            bullet_speed.animate(rate_func=rate_functions.ease_in_cubic, run_time=1.2).set_value(v/slow_down_factor),
            lag_ratio=0.8,
        )

        self.wait_until(lambda: bullet.get_x() + v/(slow_down_factor*config.frame_rate) >= pendulum_back.get_x())

        bullet.clear_updaters()
        bullet.move_to(pendulum_back)
        v_2_vector_updater(v_2_vector)

        v_2_vector.clear_updaters()
        v_3_vector = v_2_vector.copy().scale(m/(m + M)).next_to(pendulum_back, buff=0)

        self.play(AnimationGroup(
            AnimationGroup(
                ReplacementTransform(v_2_vector, v_3_vector),
                run_time=0.4, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                Flash(
                    bullet.get_right(),
                    line_length=0.1,
                    flash_radius=0.04,
                    line_stroke_width=1.5,
                    time_width=0.4,
                    color=BLACK,
                ),
                rate_func=rate_functions.ease_out_cubic,
                run_time=0.4,
            ),
        ))

        self.wait()

        circled_three = Group(
            Circle(radius=inner_radius, stroke_width=1, color=WHITE),
            MathTex('3')
        ).scale(0.5).move_to(bullet)

        bullet.set(z_index=1)
        pendulum_head.set(z_index=1)
        pendulum_back.set(z_index=1)
        self.play(AnimationGroup(
            self.camera.frame.animate(run_time=3).scale(0.9).move_to(bullet),
            circled_three.animate(rate_func=rate_functions.ease_out_cubic)
                .shift((MED_LARGE_BUFF + outer_radius - inner_radius)*DOWN),
            lag_ratio=0.3
        ))
        bullet.set(z_index=0)
        pendulum_head.set(z_index=0)
        pendulum_back.set(z_index=0)

        self.wait()

        v_2_vector.add_updater(v_2_vector_updater)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(1.5*LEFT),
            bullet.animate.shift(3*LEFT),
            ReplacementTransform(v_3_vector.copy(), v_2_vector),
            v_3_vector.animate.set_opacity(0.25),
            lag_ratio=0.05, run_time=2.5
        ))

        circled_two.set_y(circled_three.get_y())
        
        self.play(
            circled_two.animate.set_x(bullet.get_x()),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(1.5)

        v_2_vector_label.next_to(v_2_vector, buff=SMALL_BUFF)
        mass_label.move_to(bullet).set(color=BLACK)

        self.play(
            FadeIn(mass_label),
            AnimationGroup(FadeIn(v_2_vector_label, target_position=v_2_vector.get_end())),
            lag_ratio=0.5, run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        momentum_eq_2 = MathTex('p_2', '=', 'm', 'v', '_2').scale(0.4).next_to(circled_two, buff=SMALL_BUFF)
        momentum_eq_2[3:].set(color=VELOCITY_COLOR)

        mass_label_copy = mass_label.copy()

        self.play(AnimationGroup(
            FadeIn(momentum_eq_2[:2], shift=SMALL_BUFF*RIGHT),
            AnimationGroup(
                mass_label_copy.animate.move_to(momentum_eq_2[2]).set(color=WHITE),
                AnimationGroup(
                    v_2_vector_label[0].animate.move_to(momentum_eq_2[3]),
                    v_2_vector_label[1].animate.move_to(momentum_eq_2[4]),
                ),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.15, run_time=3.5
        ))

        self.add(momentum_eq_2)
        self.remove(mass_label_copy, *v_2_vector_label)

        self.wait(2)

        big_mass_label = MathTex('M', color=BLACK).scale(0.4).move_to(pendulum_back)

        v_3_vector_label = MathTex('v', '_3', color=VELOCITY_COLOR).scale(0.4).next_to(v_3_vector, RIGHT, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            v_2_vector.animate.set_opacity(0.25),
            v_3_vector.animate.set_opacity(1),
            FadeIn(big_mass_label),
            FadeIn(v_3_vector_label, target_position=v_3_vector.get_end()),
            lag_ratio=0.5, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        bullet_copy_2 = bullet.copy().set_opacity(0.25)
        self.add(bullet_copy_2)

        v_2_vector.clear_updaters()
    
        combined_mass_label = MathTex('m', '+', 'M', color=BLACK).scale(0.4).move_to(pendulum_back)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(0.75*RIGHT),
            AnimationGroup(
                bullet.animate.move_to(pendulum_back),
                mass_label.animate.move_to(combined_mass_label[0]),
            ),
            big_mass_label.animate.move_to(combined_mass_label[2]),
            FadeIn(combined_mass_label[1]),
            lag_ratio=0.05, run_time=2.5
        ))

        self.add(combined_mass_label)
        self.remove(mass_label, big_mass_label)

        self.wait()

        momentum_eq_3 = MathTex('p_3', '=', '(', 'm', '+', 'M', ')', 'v', '_3')
        momentum_eq_3.scale(0.4).next_to(circled_three, buff=SMALL_BUFF)
        momentum_eq_3[7:].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            FadeIn(momentum_eq_3[:2], shift=SMALL_BUFF*RIGHT),
            AnimationGroup(
                AnimationGroup(
                    combined_mass_label[0].animate.move_to(momentum_eq_3[3]).set(color=WHITE),
                    combined_mass_label[1].animate.move_to(momentum_eq_3[4]).set(color=WHITE),
                    combined_mass_label[2].animate.move_to(momentum_eq_3[5]).set(color=WHITE),
                ),
                AnimationGroup(
                    v_3_vector_label[0].animate.move_to(momentum_eq_3[7]),
                    v_3_vector_label[1].animate.move_to(momentum_eq_3[8]),
                ),
                AnimationGroup(
                    FadeIn(momentum_eq_3[2], shift=SMALL_BUFF*LEFT),
                    FadeIn(momentum_eq_3[6], shift=SMALL_BUFF*RIGHT),
                ),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.15, run_time=3.5
        ))

        self.add(momentum_eq_3)
        self.remove(*combined_mass_label, *v_3_vector_label)

        self.wait(2)

        max_pendulum_width = l * np.sin(delta_theta)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(
                    Group(bullet_copy_2, v_2_vector, circled_two, momentum_eq_2),
                    shift=2*LEFT,
                ),
                rate_func=rate_functions.ease_in_cubic
            ),
            self.camera.frame.animate(rate_func=rate_functions.ease_in_out_cubic)
                .scale(1.2)
                .move_to(Group(pendulum, circled_three))
                .shift(max_pendulum_width/2 * RIGHT),
            lag_ratio=0.15, run_time=3
        ))

        self.wait()

        energy_eq_3 = MathTex('E_3', '=', r'\tfrac{1}{2}', '(', 'm', '+', 'M', ')', 'v', '^2', '_3')
        energy_eq_3.scale(0.4).shift(momentum_eq_3[1].get_center() - energy_eq_3[1].get_center())
        energy_eq_3[8::2].set(color=VELOCITY_COLOR)

        momentum_eq_3_copy = momentum_eq_3.copy()
        self.remove(*momentum_eq_3)
        self.add(momentum_eq_3_copy)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(momentum_eq_3_copy[0], shift=SMALL_BUFF*UP),
                FadeIn(energy_eq_3[0], shift=SMALL_BUFF*UP),
                lag_ratio=0.05
            ),
            AnimationGroup(
                AnimationGroup(*[
                    momentum_eq_3_copy[i].animate.move_to(energy_eq_3[j])
                    for i, j in zip(range(8, 1, -1), [10, 8, 7, 6, 5, 4, 3])
                ], lag_ratio=0.05),
                FadeIn(energy_eq_3[2], shift=SMALL_BUFF*UP),
                lag_ratio=0.2
            ),
            FadeIn(energy_eq_3[9], target_position=energy_eq_3[8]),
            lag_ratio=0.1, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.add(energy_eq_3)
        self.remove(*momentum_eq_3_copy)

        self.wait()

        zero_height_line = DashedLine(
            self.camera.frame.get_left(), self.camera.frame.get_right(),
            stroke_width=1
        ).set_opacity(0.5).set_y(bullet.get_y())

        self.play(
            Create(zero_height_line),
            run_time=0.3, rate_fun=rate_functions.ease_in_cubic
        )

        self.wait(0.1)

        height_brace = BraceBetweenPoints(
            ORIGIN, 2 * (l + outer_radius) * (1 - np.cos(delta_theta)) * UP,
            buff=0.2, color=DISPLACEMENT_COLOR
        )
        height_brace.shift(pendulum_back.get_right())

        height_brace.scale(0.5, about_point=height_brace.get_corner(DL))
        height_brace.stretch_to_fit_height(1e-6, about_point=height_brace.get_bottom())

        self.add(height_brace)

        pendulum_copy_3 = pendulum.copy()
        for mob in pendulum_copy_3.submobjects:
            mob.set_opacity(0.25)
        bullet_copy_3 = bullet.copy().set_opacity(0.25)
        v_3_vector_copy = v_3_vector.copy().set_opacity(0.25)

        self.add(pendulum_copy_3, bullet_copy_3, v_3_vector_copy)

        height_brace.add_updater(
            lambda mob: mob.set_x(pendulum_back.get_right()[0] + 0.2)
                .stretch_to_fit_height(
                    max(bullet.get_y() - bullet_copy_3.get_y(), 1e-6),
                    about_point=height_brace.get_bottom()
                )
        )

        to_v_3_start = v_3_vector.get_start() - pendulum_shaft.get_start()
        to_v_3_end = v_3_vector.get_end() - pendulum_shaft.get_start() + 0.01*RIGHT
        v_3_angle_diff = angle_between_vectors(to_v_3_start, to_v_3_end)

        self.play(
            Rotate(
                Group(pendulum, bullet),
                delta_theta,
                about_point=pendulum_shaft.get_start()
            ),
            Rotate(
                v_3_vector,
                delta_theta - v_3_angle_diff,
                about_point=pendulum_shaft.get_start()
            ),
            run_time=1.5, rate_func=rate_functions.ease_in_out_sine
        )

        height_brace.clear_updaters()
        self.remove(v_3_vector)

        h_label = MathTex('h', color=DISPLACEMENT_COLOR).scale(0.4).next_to(height_brace, buff=SMALL_BUFF)

        self.play(
            FadeIn(h_label, target_position=height_brace),
            run_time=0.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        circled_four = Group(
            Circle(radius=inner_radius, stroke_width=1, color=WHITE),
            MathTex('4')
        ).scale(0.5).move_to(bullet)

        bullet.set(z_index=1)
        pendulum_head.set(z_index=1)
        pendulum_back.set(z_index=0)
        self.play(AnimationGroup(
            circled_four.animate.shift((MED_LARGE_BUFF + outer_radius - inner_radius)*UP),
            rate_func=rate_functions.ease_out_cubic
        ))
        bullet.set(z_index=0)
        pendulum_head.set(z_index=0)
        pendulum_back.set(z_index=0)

        energy_eq_4 = MathTex('E_4', '=', '(', 'm', '+', 'M', ')', 'g', 'h')
        energy_eq_4.scale(0.4).next_to(circled_four, buff=SMALL_BUFF)
        energy_eq_4[8].set(color=DISPLACEMENT_COLOR)

        self.play(AnimationGroup(
            FadeIn(energy_eq_4[:2], shift=SMALL_BUFF*RIGHT),
            FadeIn(energy_eq_4[3:6], target_position=bullet),
            AnimationGroup(
                FadeIn(energy_eq_4[2], shift=SMALL_BUFF*LEFT),
                FadeIn(energy_eq_4[6], shift=SMALL_BUFF*RIGHT),
            ),
            FadeIn(energy_eq_4[7], shift=SMALL_BUFF*DOWN),
            AnimationGroup(
                h_label.animate.move_to(energy_eq_4[8]),
                FadeOut(zero_height_line),
                lag_ratio=0.05
            ),
            lag_ratio=0.1, run_time=5.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.add(energy_eq_4)
        self.remove(h_label)

        circled_num_buff = circled_four.get_x() - circled_three.get_x()
        spring_shift_x_val = circled_three.get_x() - 2*circled_num_buff - circled_one.get_x()
        Group(spring, spring_brace, circled_one).shift(spring_shift_x_val*RIGHT)

        bullet_copy_2.set_x(mid(circled_one.get_x(), circled_three.get_x()))
        v_2_vector.next_to(bullet_copy_2, buff=0)
        circled_two.next_to(bullet_copy_2, UP, buff=SMALL_BUFF)

        energy_eq_3_copy = energy_eq_3.copy()
        energy_eq_4_copy = energy_eq_4.copy()

        energy_eq_1.scale(1.2).next_to(circled_one, DOWN, buff=1.8)
        energy_eq_2.scale(1.2).set_x(circled_two.get_x()).set_y(energy_eq_1.get_y())
        energy_eq_3.scale(1.2).set_x(circled_three.get_x()).set_y(energy_eq_1.get_y()),
        energy_eq_4.scale(1.2).set_x(circled_four.get_x()).set_y(energy_eq_1.get_y()),

        momentum_eq_2.scale(1.2).next_to(energy_eq_2, DOWN)
        momentum_eq_3.scale(1.2).set_x(energy_eq_3.get_x()).set_y(momentum_eq_2.get_y())

        self.add(momentum_eq_2, momentum_eq_3)
        self.remove(*energy_eq_3, *energy_eq_4)
        self.add(energy_eq_3_copy, energy_eq_4_copy)

        self.play(AnimationGroup(
            self.camera.frame.animate.set(height=8).scale(0.85).move_to(Group(
                spring, pendulum, spring_brace, height_brace,
                circled_one, circled_two, circled_three, circled_four,
                energy_eq_1, energy_eq_2, energy_eq_3, energy_eq_4,
                momentum_eq_2, momentum_eq_3
            )),
            FadeIn(
                Group(bullet_copy_2, v_2_vector, circled_two),
                shift=RIGHT,
                rate_func=rate_functions.ease_out_cubic
            ),
            ReplacementTransform(
                energy_eq_3_copy, energy_eq_3,
                rate_func=rate_functions.ease_in_out_cubic
            ),
            ReplacementTransform(
                energy_eq_4_copy, energy_eq_4,
                rate_func=rate_functions.ease_in_out_cubic
            ),
            lag_ratio=0.05, run_time=5
        ))

        self.wait(2)

        eq_sign_1_2 = MathTex('=', color=EQUALITY_COLOR).scale(0.6).set_y(energy_eq_1.get_y())
        eq_sign_2_3 = eq_sign_1_2.copy().set_y(momentum_eq_2.get_y())
        eq_sign_3_4 = eq_sign_1_2.copy()

        eq_sign_1_2.set_x(mid(energy_eq_1.get_right()[0], energy_eq_2.get_left()[0]))
        eq_sign_2_3.set_x(mid(momentum_eq_2.get_right()[0], momentum_eq_3.get_left()[0]))
        eq_sign_3_4.set_x(mid(energy_eq_3.get_right()[0], energy_eq_4.get_left()[0]))

        not_eq_sign_2_3 = MathTex(r'\neq', color=EQUALITY_COLOR).scale(0.6).set_y(energy_eq_1.get_y())
        not_eq_sign_2_3.set_x(mid(energy_eq_2.get_right()[0], energy_eq_3.get_left()[0]))

        self.play(AnimationGroup(
            Write(eq_sign_1_2),
            Write(eq_sign_2_3),
            Write(eq_sign_3_4),
            lag_ratio=0.5, run_time=2
        ))

        self.wait(1.5)

        self.play(Write(not_eq_sign_2_3), run_time=1.2)
        
        self.wait()
