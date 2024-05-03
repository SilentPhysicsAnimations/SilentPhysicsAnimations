from manim import *

#####################

CHARGE_COLOR = YELLOW
VELOCITY_COLOR = BLUE

#####################

class VelocitySelectorScene(MovingCameraScene):
    def construct(self):
        positron = Dot(radius=0.1, color=CHARGE_COLOR, z_index=100).shift(10*LEFT)

        self.add(positron)

        fw, fh = self.camera.frame_width, self.camera.frame_height

        fill = { 'color': '#ababab', 'stroke_width': 0, 'fill_opacity': 1 }
        top = Rectangle(width=18/16*fw, height=fh/5, **fill).to_edge(UP, buff=0)
        bottom = top.copy().to_edge(DOWN, buff=0)
        back = Rectangle(width=top.width, height=fh - top.height - bottom.height, **fill)
        back.set(color='#242424')

        top.set_sheen(-0.5, UP)
        bottom.set_sheen(-0.5, DOWN)
        back.set_sheen(-0.25, UP)

        self.add(top, bottom, back)

        m, q = 1, 1
        E, B = 8, 1
        v = (E/B)*RIGHT

        elec_length = 1
        vector_config = {
            'max_tip_length_to_length_ratio': 0.15,
            'max_stroke_width_to_length_ratio': 100,
        }

        velocity_vec = Vector(**vector_config, color=VELOCITY_COLOR, z_index=1).set_opacity(0.95)
        velocity_vec.add_updater(lambda mob: mob.put_start_and_end_on(ORIGIN, v * elec_length/E).shift(positron.get_center()))

        self.add(velocity_vec)

        self.play(
            positron.animate.center(),
            run_time=4, rate_func=rate_functions.ease_out_sine
        )

        self.wait(2)

        marks = np.linspace(-fw/2, fw/2, 16)

        pluses = Group(*[MathTex('+', color=BLACK).scale(0.8).set_x(x) for x in marks]).shift((top.get_bottom()[1] + 0.33) * UP)
        minuses = Group(*[MathTex('-', color=BLACK).scale(0.8).set_x(x) for x in marks]).shift((bottom.get_top()[1] - 0.33) * UP)

        electrical_field = Group(*[Vector((back.height - 0.66)*DOWN, stroke_width=10).move_to(x*RIGHT).set_opacity(0.2) for x in marks])

        self.play(AnimationGroup(
            *[
                AnimationGroup(
                    FadeIn(plus, shift=0.25*DOWN),
                    FadeIn(minus, shift=0.25*UP),
                    GrowArrow(vec),
                    rate_func=rate_functions.ease_out_cubic
                )
                for plus, minus, vec in zip(pluses, minuses, electrical_field)
            ],
            lag_ratio=0.01, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        electric_force_vec = Vector(elec_length*DOWN, **vector_config).set_opacity(0.95)
        electric_force_eq = MathTex(r'\vec{F_E}', '=', 'q', r'\vec{E}', background_stroke_width=4)
        electric_force_eq.scale(0.7).next_to(electric_force_vec, buff=0.1)
        electric_force_eq[2].set(color=CHARGE_COLOR)

        self.play(AnimationGroup(
            GrowArrow(electric_force_vec),
            FadeIn(electric_force_eq[0], shift=0.1*RIGHT),
            lag_ratio=0.2, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.2)

        self.play(AnimationGroup(
            FadeIn(electric_force_eq[1]),
            FadeIn(electric_force_eq[2], target_position=positron),
            FadeIn(electric_force_eq[3]),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        self.play(AnimationGroup(
            FadeOut(electric_force_eq),
            rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        time_speed = 0.2

        B = 0

        F = np.zeros(3)
        def positron_updater(mob, dt):
            nonlocal v, F
            dt *= time_speed
            F = q * (np.cross(v, B * IN) + E * DOWN)
            v += F/m * dt
            mob.shift(v * dt)

        positron.add_updater(positron_updater)

        electric_force_vec.add_updater(lambda mob: mob.put_start_and_end_on(ORIGIN, elec_length*DOWN).shift(positron.get_center()))

        self.wait_until(lambda: positron.get_bottom()[1] <= bottom.get_top()[1])

        positron.clear_updaters()
        velocity_vec.clear_updaters()
        electric_force_vec.clear_updaters()

        self.wait(3)

        self.play(AnimationGroup(
            FadeOut(electrical_field),
            FadeOut(pluses),
            FadeOut(minuses),
            FadeOut(electric_force_vec),
            positron.animate.move_to(DOWN),
            velocity_vec.animate.put_start_and_end_on(ORIGIN, RIGHT).shift(DOWN),
            run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        crossed_circle = VGroup(
            Circle(radius=0.2, color=WHITE, stroke_opacity=0.2),
            Line(0.18*LEFT, 0.18*RIGHT, stroke_opacity=0.2).rotate(45*DEGREES),
            Line(0.18*LEFT, 0.18*RIGHT, stroke_opacity=0.2).rotate(-45*DEGREES)
        )
        magnetic_field = Group()
        for x in marks + (marks[1] - marks[0])/2:
            for y in np.linspace(bottom.get_top()[1] + 0.89, top.get_bottom()[1] - 0.89, 4):
                magnetic_field.add(crossed_circle.copy().move_to([x, y, 0]))

        self.play(
            FadeIn(magnetic_field),
            lag_ratio=0.0025, run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.5)

        magnetic_force_vec = Vector(UP, **vector_config).set_opacity(0.95).shift(DOWN)
        magnetic_force_eq = MathTex(r'\vec{F_B}', '=', 'q', r'\vec{v}', r'\times', r'\vec{B}', background_stroke_width=4)
        magnetic_force_eq.scale(0.7).next_to(magnetic_force_vec, buff=0.1)
        magnetic_force_eq[2].set(color=CHARGE_COLOR)
        magnetic_force_eq[3].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            GrowArrow(magnetic_force_vec),
            FadeIn(magnetic_force_eq[0], shift=0.1*RIGHT),
            lag_ratio=0.2, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.2)

        self.play(AnimationGroup(
            FadeIn(magnetic_force_eq[1]),
            FadeIn(magnetic_force_eq[2], target_position=positron),
            FadeIn(magnetic_force_eq[3], target_position=velocity_vec),
            FadeIn(magnetic_force_eq[4]),
            FadeIn(magnetic_force_eq[5]),
            lag_ratio=0.1, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        self.play(AnimationGroup(
            FadeOut(magnetic_force_eq),
            rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        self.play(Rotate(
            Group(positron, velocity_vec, magnetic_force_vec),
            angle=2*PI, about_point=ORIGIN,
            run_time=2*PI/1.4, rate_func=rate_functions.linear
        ))

        self.wait(2.5)

        electric_force_vec.put_start_and_end_on(ORIGIN, 1e-6*DOWN).move_to(positron)
        self.add(electric_force_vec)

        pos_pos = (marks[2] - marks[0])*LEFT

        self.play(AnimationGroup(
            AnimationGroup(*[
                AnimationGroup(
                    FadeIn(plus, shift=0.25*DOWN),
                    FadeIn(minus, shift=0.25*UP),
                    GrowArrow(vec),
                    rate_func=rate_functions.ease_out_cubic
                )
                for plus, minus, vec in zip(pluses, minuses, electrical_field)
            ], lag_ratio=0.01, rate_func=rate_functions.ease_out_cubic),
            AnimationGroup(
                positron.animate.move_to(pos_pos),
                electric_force_vec.animate.put_start_and_end_on(ORIGIN, elec_length*DOWN).shift(pos_pos),
                magnetic_force_vec.animate.put_start_and_end_on(ORIGIN, elec_length*UP).shift(pos_pos),
                velocity_vec.animate.put_start_and_end_on(ORIGIN, elec_length*RIGHT).shift(pos_pos),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.1, run_time=4
        ))

        self.wait(2)

        electric_force_eq.next_to(electric_force_vec, buff=0.1)
        magnetic_force_eq.next_to(magnetic_force_vec, buff=0.1)

        self.play(AnimationGroup(
            FadeIn(magnetic_force_eq, shift=0.1*RIGHT),
            FadeIn(electric_force_eq, shift=0.1*RIGHT),
            lag_ratio=0.2, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v_eq_0 = MathTex(r'\vec{F_B}', '+', r'\vec{F_E}', '=', r'\vec{0}', background_stroke_width=4)
        v_eq_0.scale(0.7).next_to(velocity_vec)

        self.play(AnimationGroup(
            FadeIn(v_eq_0[0], target_position=magnetic_force_eq[0]),
            FadeIn(v_eq_0[1], shift=0.1*RIGHT),
            FadeIn(v_eq_0[2], target_position=electric_force_eq[0]),
            FadeIn(v_eq_0[3], shift=0.1*RIGHT),
            FadeIn(v_eq_0[4], shift=0.1*RIGHT),
            lag_ratio=0.05, run_time=4.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v_eq_1 = MathTex('q', r'\vec{v}', r'\times', r'\vec{B}', '+', 'q', r'\vec{E}', '=', r'\vec{0}', background_stroke_width=4)
        v_eq_1.scale(0.7).next_to(velocity_vec)
        v_eq_1[0:6:5].set(color=CHARGE_COLOR)
        v_eq_1[1].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            FadeOut(v_eq_0[0]),
            FadeOut(v_eq_0[2]),
            FadeOut(magnetic_force_eq[:2]),
            ReplacementTransform(magnetic_force_eq[2:], v_eq_1[:4]),
            ReplacementTransform(v_eq_0[1], v_eq_1[4]),
            FadeOut(electric_force_eq[:2]),
            ReplacementTransform(electric_force_eq[2:], v_eq_1[5:7]),
            ReplacementTransform(v_eq_0[3], v_eq_1[7]),
            ReplacementTransform(v_eq_0[4], v_eq_1[8]),
            lag_ratio=0.05, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v_eq_2 = MathTex('q', 'v', 'B', '-', 'q', 'E', '=', '0', background_stroke_width=4)
        v_eq_2.scale(0.7).next_to(velocity_vec)
        v_eq_2[0:5:4].set(color=CHARGE_COLOR)
        v_eq_2[1].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            TransformMatchingShapes(v_eq_1[:4], v_eq_2[:3], fade_transform_mismatches=True),
            ReplacementTransform(v_eq_1[4], v_eq_2[3]),
            TransformMatchingShapes(v_eq_1[5:7], v_eq_2[4:6], fade_transform_mismatches=True),
            ReplacementTransform(v_eq_1[7], v_eq_2[6]),
            TransformMatchingShapes(v_eq_1[8], v_eq_2[7], fade_transform_mismatches=True),
            lag_ratio=0.05, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v_eq_3 = MathTex('q', 'v', 'B', '=', 'q', 'E', background_stroke_width=4)
        v_eq_3.scale(0.7).next_to(velocity_vec)
        v_eq_3[0:5:4].set(color=CHARGE_COLOR)
        v_eq_3[1].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(v_eq_2[:3], v_eq_3[:3]),
            FadeOut(v_eq_2[3]),
            ReplacementTransform(v_eq_2[4:6], v_eq_3[4:6]),
            ReplacementTransform(v_eq_2[6], v_eq_3[3]),
            FadeOut(v_eq_2[7]),
            lag_ratio=0.05, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        self.play(AnimationGroup(
            FadeOut(v_eq_3[0:5:4], shift=0.1*DOWN),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.2)

        v_eq = MathTex('v', '=', '{E', r'\over', 'B}', background_stroke_width=4)
        v_eq.scale(0.7).next_to(velocity_vec)
        v_eq[0].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(v_eq_3[1], v_eq[0]),
            ReplacementTransform(v_eq_3[3], v_eq[1]),
            ReplacementTransform(v_eq_3[5], v_eq[2]),
            FadeIn(v_eq[3], shift=0.1*UP),
            ReplacementTransform(v_eq_3[2], v_eq[4]),
            lag_ratio=0.05, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(3)

        positron_copy = positron.copy()
        positron_middle_path = TracedPath(
            positron_copy.get_center,
            stroke_color=CHARGE_COLOR, stroke_opacity=0.5, stroke_width=6
        )

        self.add(positron_copy, positron_middle_path)

        self.play(
            positron_copy.animate.shift(RIGHT),
            rate_func=rate_functions.linear
        )

        self.play(
            Group(positron_copy, v_eq).animate.shift((2*abs(pos_pos[0]) - 1)*RIGHT),
            run_time=2*abs(pos_pos[0])-1, rate_func=rate_functions.linear
        )

        positron_middle_path.clear_updaters()
        positron_middle_path.become(
            Line(pos_pos, -pos_pos, stroke_color=CHARGE_COLOR, stroke_opacity=0.5, stroke_width=6)
        )

        self.wait()

        self.play(AnimationGroup(
            FadeOut(positron_copy),
            v_eq.animate.set_opacity(0.5),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            velocity_vec.animate.scale(1.5, about_point=positron.get_center()),
            magnetic_force_vec.animate.scale(1.5, about_point=positron.get_center()),
            run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        positron_copy.move_to(positron)
        positron_upper_path = TracedPath(
            positron_copy.get_center,
            stroke_color=CHARGE_COLOR, stroke_opacity=0.5, stroke_width=6
        )

        self.add(positron_copy, positron_upper_path)

        E, B = 1.5, 1
        time_speed = B/E

        v = 1.5*(E/B)*RIGHT
        F = np.zeros(3)
        positron_copy.add_updater(positron_updater)

        self.wait_until(lambda: positron_copy.get_x() >= abs(pos_pos[0]))

        positron_copy.clear_updaters()
        positron_upper_path.clear_updaters()

        self.wait()

        v_bigger = MathTex('v', '>', '{E', r'\over', 'B}', background_stroke_width=4)
        v_bigger.scale(0.7).next_to(positron_copy, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER-0.1)
        v_bigger[0].set(color=VELOCITY_COLOR)

        self.play(
            TransformMatchingTex(v_eq.copy(), v_bigger),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(AnimationGroup(
            FadeOut(positron_copy),
            v_bigger.animate.set_opacity(0.5),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            velocity_vec.animate.scale(1/3, about_point=positron.get_center()),
            magnetic_force_vec.animate.scale(1/3, about_point=positron.get_center()),
            run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(0.5)

        positron_copy.move_to(positron)
        positron_lower_path = TracedPath(
            positron_copy.get_center,
            stroke_color=CHARGE_COLOR, stroke_opacity=0.5, stroke_width=6
        )

        self.add(positron_copy, positron_lower_path)

        v = 0.5*(E/B)*RIGHT
        F = np.zeros(3)
        positron_copy.add_updater(positron_updater)

        self.wait_until(lambda: positron_copy.get_x() >= abs(pos_pos[0]))

        positron_copy.clear_updaters()
        positron_lower_path.clear_updaters()

        self.wait()

        v_smaller = MathTex('v', '<', '{E', r'\over', 'B}', background_stroke_width=4)
        v_smaller.scale(0.7).next_to(positron_copy, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER - 0.1)
        v_smaller[0].set(color=VELOCITY_COLOR)

        self.play(
            TransformMatchingTex(v_eq.copy(), v_smaller),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(AnimationGroup(
            FadeOut(positron_copy),
            v_smaller.animate.set_opacity(0.5),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        wall_top = Rectangle(width=marks[1]-marks[0], height=(back.height-0.06)/2, **fill).next_to(top, DOWN, buff=0)
        wall_bottom = wall_top.copy().next_to(bottom, UP, buff=0)
        Group(wall_top, wall_bottom).move_to((marks[-1] + marks[1] - marks[0])*RIGHT, aligned_edge=LEFT).shift(0.06*LEFT)

        arc_top = Sector(**fill).stretch_to_fit_width(wall_top.width).stretch_to_fit_height(top.height).next_to(wall_top, UP, buff=0)
        arc_bottom = arc_top.copy().flip(axis=RIGHT).next_to(wall_bottom, DOWN, buff=0)

        arc_top.set_sheen(-0.5, UP)
        arc_bottom.set_sheen(-0.5, DOWN)

        self.add(wall_top, wall_bottom, arc_top, arc_bottom)

        self.play(AnimationGroup(
            v_eq.animate.set_opacity(1).next_to(wall_bottom).align_to(wall_bottom, UP).shift(DEFAULT_MOBJECT_TO_MOBJECT_BUFFER*DOWN),
            positron_middle_path.animate.put_start_and_end_on(positron.get_center(), 18*RIGHT),
            self.camera.frame.animate.move_to(Group(wall_top, wall_bottom).get_left()),
            AnimationGroup(
                FadeOut(positron_upper_path),
                FadeOut(positron_lower_path),
                FadeOut(v_bigger),
                FadeOut(v_smaller),
                FadeOut(velocity_vec),
                FadeOut(electric_force_vec),
                FadeOut(magnetic_force_vec),
            ),
            lag_ratio=0.05, run_time=12, rate_func=rate_functions.ease_in_out_sine
        ))

        self.wait()

        E, B = 4, 1
        time_speed = 2

        def positron_updater(v_0, will_hit_wall):
            v = v_0
            F = np.zeros(3)
            def updater(mob, dt):
                nonlocal v, F
                dt *= time_speed
                F = q * (np.cross(v, B * IN) + E * DOWN)
                v += F/m * dt
                mob.shift(v * dt)
                if will_hit_wall and (mob.get_x() >= wall_top.get_left()[0] or not bottom.get_top()[1] < mob.get_y() < top.get_bottom()[1]):
                    mob.clear_updaters()
            return updater

        for _ in range(200):
            positron = Dot(radius=0.03, color=CHARGE_COLOR)
            positron_trail = TracedPath(
                positron.get_center, stroke_color=CHARGE_COLOR,
                dissipating_time=0.2, stroke_opacity=[1, 0]
            )

            hits_wall = True
            v_0 = np.random.uniform(0.5, 2)
            if abs(v_0 - 1) < 0.1:
                v_0 = 1
                hits_wall = False
            v_0 *= (E/B) * RIGHT
            positron.add_updater(positron_updater(v_0, hits_wall))

            self.add(positron, positron_trail)
            self.wait(0.05)

