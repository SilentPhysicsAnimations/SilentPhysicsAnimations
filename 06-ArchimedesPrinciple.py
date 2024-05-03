from manim import *

#######################################

WATER_DENSITY_LABEL = r'\rho_{\text{w}}'
WATER_DENSITY_COLOR = BLUE

OBJECT_DENSITY_LABEL = r'\rho_{\text{obj}}'
OBJECT_DENSITY_COLOR = WHITE

GRAVITY_FORCE_LABEL = 'F_g'
GRAVITY_FORCE_COLOR = GREEN

BUOYANT_FORCE_LABEL = 'F_b'
BUOYANT_FORCE_COLOR = RED

VOLUME_COLOR = YELLOW

#######################################

TEX = { 'tex_template': TexTemplateLibrary.default }

#######################################

class ArchimedesPrincipleScene(MovingCameraScene):
    def construct(self):
        R = 1
        r = 1/3
        h = 1.5*r

        glass_buff = 0.05
        water_width = 2*R
        height_delta = (r/R)**2 * h
        
        water = Rectangle(
            width=2*R, height=3*R,
            stroke_width=0, fill_color=BLUE, fill_opacity=0.5,
            sheen_factor=0.2, sheen_direction=LEFT,
            z_index=1
        )
        glass = VMobject(stroke_width=2)
        glass.start_new_path(water.get_corner(UL) + glass_buff*LEFT + height_delta*UP)
        glass.add_line_to(water.get_corner(DL) + glass_buff*DL)
        glass.add_line_to(water.get_corner(DR) + glass_buff*DR)
        glass.add_line_to(water.get_corner(UR) + glass_buff*RIGHT + height_delta*UP)

        Group(glass, water).center()

        self.add(glass, water)

        obj = Rectangle(
            width=2*r, height=h,
            stroke_width=0, fill_opacity=1,
            sheen_factor=-0.3, sheen_direction=RIGHT
        ).next_to(water, UP)

        self.add(obj)

        water_raise_bottom_bounds = (water.get_top()[1], glass.get_top()[1] - h)
        def water_updater(mob):
            ratio = (obj.get_bottom()[1] - water_raise_bottom_bounds[0]) / (water_raise_bottom_bounds[1] - water_raise_bottom_bounds[0])
            mob.stretch_to_fit_height(3*R + height_delta * np.clip(ratio, 0, 1), about_edge=DOWN)

        water.add_updater(water_updater)

        self.play(obj.animate.center(), run_time=4)

        water.clear_updaters()

        self.wait(2)

        self.play(
            self.camera.frame.animate.scale(0.67),
            water.animate.set_opacity(0.25),
        )

        self.wait()

        h_1_line = DashedLine(water.get_top(), obj.get_top(), color=VOLUME_COLOR, stroke_width=1.5).shift((r - 0.0075)*LEFT)
        h_2_line = DashedLine(water.get_top(), obj.get_bottom(), color=VOLUME_COLOR, stroke_width=1.5).shift((r + 0.0375)*RIGHT)
        h_1_label = MathTex('h_1', color=VOLUME_COLOR).scale(0.5).next_to(h_1_line, LEFT, buff=SMALL_BUFF)
        h_2_label = MathTex('h_2', color=VOLUME_COLOR).scale(0.5).next_to(h_2_line, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            AnimationGroup(
                Create(h_1_line),
                FadeIn(h_1_label, target_position=h_1_line),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                Create(h_2_line),
                FadeIn(h_2_label, target_position=h_2_line),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.7, run_time=3
        ))

        self.wait(2)

        pressure_vec = Vector(
            0.25*DOWN, color=BUOYANT_FORCE_COLOR,
            max_tip_length_to_length_ratio=0.4,
            max_stroke_width_to_length_ratio=10
        )

        vec_buff = (obj.width - 5*pressure_vec.width) / 4

        P_1_vecs = VGroup(*[pressure_vec.copy() for _ in range(0, 2)])
        P_1_vecs.arrange(buff=2*vec_buff+pressure_vec.width).next_to(obj, UP, buff=0.05)

        P_1_eq = MathTex('P_1', '=', WATER_DENSITY_LABEL, 'g', 'h_1').scale(0.5)
        P_1_eq[0].set(color=BUOYANT_FORCE_COLOR)
        P_1_eq[2].set(color=WATER_DENSITY_COLOR)
        P_1_eq[4].set(color=VOLUME_COLOR)
        P_1_eq.next_to(glass).set_y(P_1_vecs.get_y())

        self.play(AnimationGroup(
            AnimationGroup(
                *[FadeIn(vec, shift=0.25*DOWN) for vec in P_1_vecs],
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(P_1_eq[:4], shift=0.2*RIGHT),
                h_1_label.animate.move_to(P_1_eq[4]),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            run_time=4
        ))

        self.remove(h_1_label)
        self.add(P_1_eq)

        self.wait(0.3)

        P_2_vecs = VGroup(*[pressure_vec.copy().rotate(PI) for _ in range(0, 5)])
        P_2_vecs.arrange(buff=vec_buff).next_to(obj, DOWN, buff=0.05)

        P_2_eq = MathTex('P_2', '=', WATER_DENSITY_LABEL, 'g', 'h_2').scale(0.5)
        P_2_eq[0].set(color=BUOYANT_FORCE_COLOR)
        P_2_eq[2].set(color=WATER_DENSITY_COLOR)
        P_2_eq[4].set(color=VOLUME_COLOR)
        P_2_eq.next_to(glass).set_y(P_2_vecs.get_y())

        self.play(AnimationGroup(
            AnimationGroup(
                *[FadeIn(vec, shift=0.25*UP) for vec in P_2_vecs],
                lag_ratio=0.067, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(P_2_eq[:4], shift=0.2*RIGHT),
                h_2_label.animate.move_to(P_2_eq[4]),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            run_time=4
        ))

        self.remove(h_2_label)
        self.add(P_2_eq)

        self.wait(2)

        P_diff_eq = MathTex(
            'P', '&=', 'P_2', '-', 'P_1',
            r'\\&=', WATER_DENSITY_LABEL, 'g', 'h_2', '-', WATER_DENSITY_LABEL, 'g', 'h_1'
        ).scale(0.5).next_to(glass)
        P_diff_eq[0:5:2].set(color=BUOYANT_FORCE_COLOR)
        P_diff_eq[6::4].set(color=WATER_DENSITY_COLOR)
        P_diff_eq[8::4].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(RIGHT),
            AnimationGroup(
                FadeOut(P_1_vecs, shift=0.05*DOWN),
                FadeOut(P_2_vecs[1::2], shift=0.05*UP),
            ),
            AnimationGroup(
                P_2_vecs[0].animate.shift(pressure_vec.width/2*RIGHT),
                P_2_vecs[4].animate.shift(pressure_vec.width/2*LEFT),
            ),
            P_1_eq[1].animate.move_to(P_diff_eq[1]),
            P_2_eq[0].animate.move_to(P_diff_eq[2]),
            P_1_eq[0].animate.move_to(P_diff_eq[4]),
            P_2_eq[1].animate.move_to(P_diff_eq[5]),
            P_2_eq[2:].animate.move_to(P_diff_eq[6:9]),
            P_1_eq[2:].animate.move_to(P_diff_eq[10:]),
            FadeIn(P_diff_eq[0]),
            FadeIn(P_diff_eq[3]),
            FadeIn(P_diff_eq[9]),
            lag_ratio=0.05, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*P_1_eq, *P_2_eq)
        self.add(P_diff_eq)

        self.wait()

        P_diff_eq_2 = MathTex(
            'P', '&=', 'P_2', '-', 'P_1',
            r'\\&=', WATER_DENSITY_LABEL, 'g', '(', 'h_2', '-', 'h_1', ')'
        ).scale(0.5).next_to(glass)
        P_diff_eq_2[0:5:2].set(color=BUOYANT_FORCE_COLOR)
        P_diff_eq_2[6].set(color=WATER_DENSITY_COLOR)
        P_diff_eq_2[9::2].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            P_diff_eq[10:12].animate.move_to(P_diff_eq_2[6:8]),
            P_diff_eq[8].animate.move_to(P_diff_eq_2[9]),
            P_diff_eq[9].animate.move_to(P_diff_eq_2[10]),
            P_diff_eq[12].animate.move_to(P_diff_eq_2[11]),
            FadeIn(P_diff_eq_2[8::4]),
            lag_ratio=0.05, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*P_diff_eq)
        self.add(P_diff_eq_2)

        self.wait(0.8)

        line_top = ValueTracker(water.get_top()[1])

        def line_updater(mob):
            for dash in mob.submobjects:
                start, end, y = dash.get_start(), dash.get_end(), line_top.get_value()
                if end[1] < y < start[1]:
                    start[1] = y
                    dash.put_start_and_end_on(start, end)
                elif y <= end[1]:
                    dash.set(stroke_width=0)

        h_1_line.add_updater(line_updater)
        h_2_line.add_updater(line_updater)

        self.play(AnimationGroup(
            line_top.animate.set_value(obj.get_top()[1]-0.05),
            run_time=2, rate_func=rate_functions.ease_out_quad
        ))

        h_1_line.clear_updaters()
        h_2_line.clear_updaters()
        self.remove(*h_1_line)

        h_line = DashedLine(obj.get_top(), obj.get_bottom(), color=VOLUME_COLOR, stroke_width=1.5).shift((r + 0.0375)*RIGHT)
        h_label = MathTex('h', color=VOLUME_COLOR).scale(0.5).next_to(h_line, RIGHT, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            ReplacementTransform(h_2_line[-len(h_line):], h_line),
            FadeIn(h_label, shift=SMALL_BUFF*RIGHT),
            lag_ratio=0.1, run_time=1.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*h_2_line)

        self.wait(1.5)

        P_diff_eq_3 = MathTex(
            'P', '&=', 'P_2', '-', 'P_1',
            r'\\&=', WATER_DENSITY_LABEL, 'g', 'h',
        ).scale(0.5).next_to(glass)
        P_diff_eq_3[0:5:2].set(color=BUOYANT_FORCE_COLOR)
        P_diff_eq_3[6].set(color=WATER_DENSITY_COLOR)
        P_diff_eq_3[8].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            FadeOut(P_diff_eq_2[8:], shift=0.25*RIGHT),
            h_label.animate.move_to(P_diff_eq_3[8]),
            FadeOut(h_line),
            lag_ratio=0.05, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*P_diff_eq_2, h_label)
        self.add(P_diff_eq_3)

        self.wait()

        P_diff_eq_4 = MathTex(
            'P', '=', WATER_DENSITY_LABEL, 'g', 'h'
        ).scale(0.5).next_to(glass)
        P_diff_eq_4[0].set(color=BUOYANT_FORCE_COLOR)
        P_diff_eq_4[2].set(color=WATER_DENSITY_COLOR)
        P_diff_eq_4[4].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(LEFT),
            FadeOut(P_diff_eq_3[1:5], shift=SMALL_BUFF*UP),
            P_diff_eq_3[0].animate.move_to(P_diff_eq_4[0]),
            P_diff_eq_3[5:].animate.move_to(P_diff_eq_4[1:]),
            lag_ratio=0.05, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*P_diff_eq_3)
        self.add(P_diff_eq_4)

        self.wait(2)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(LEFT),
            P_2_vecs[0::2].animate.set_opacity(0.5),
            P_diff_eq_4.animate.set_opacity(0.5),
            lag_ratio=0.05, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        bottom = Circle(radius=obj.width/2, fill_opacity=1, color=VOLUME_COLOR)
        bottom.stretch_to_fit_height(1e-6).move_to(obj.get_bottom())

        self.play(
            Create(bottom),
            run_time=1.5, rate_func=rate_functions.linear
        )

        self.wait(1.5)

        self.play(
            bottom.animate.stretch_to_fit_height(obj.width).next_to(glass, LEFT),
            run_time=2.5, rate_func=rate_functions.ease_in_out_sine
        )

        self.wait()

        A_label = MathTex('A', color=VOLUME_COLOR).scale(0.5).move_to(bottom)
        bottom.rotate(90*DEGREES)

        self.play(
            ReplacementTransform(bottom, A_label),
            run_time=1.5,
        )

        self.wait(2)

        force_eq = MathTex(BUOYANT_FORCE_LABEL, '=', 'P', 'A').scale(0.5)
        force_eq[0].set(color=BUOYANT_FORCE_COLOR)
        force_eq[2].set(color=BUOYANT_FORCE_COLOR)
        force_eq[3].set(color=VOLUME_COLOR)
        force_eq.next_to(P_diff_eq_4, UP, aligned_edge=LEFT, buff=0.15)

        force_vec = Vector(
            0.6*UP, color=BUOYANT_FORCE_COLOR,
            max_tip_length_to_length_ratio=0.3,
            max_stroke_width_to_length_ratio=50,
        ).next_to(obj, DOWN, buff=0.05)

        self.play(AnimationGroup(
            self.camera.frame.animate.shift(RIGHT),
            AnimationGroup(
                AnimationGroup(*[Transform(vec, force_vec) for vec in P_2_vecs[0::2]]),
                FadeIn(force_eq[0], shift=SMALL_BUFF*RIGHT),
                lag_ratio=0.1
            ),
            AnimationGroup(
                FadeIn(force_eq[1], shift=SMALL_BUFF*RIGHT),
                ReplacementTransform(P_diff_eq_4[0].copy().set_opacity(0), force_eq[2]),
                ReplacementTransform(A_label, force_eq[3]),
                lag_tatio=0.1,
            ),
            lag_ratio=0.67, run_time=5
        ))

        self.remove(*P_2_vecs)
        self.add(force_vec)

        self.wait(2)

        force_eq_2 = MathTex(BUOYANT_FORCE_LABEL, '=', WATER_DENSITY_LABEL, 'g', 'h', 'A').scale(0.5)
        force_eq_2[0].set(color=BUOYANT_FORCE_COLOR)
        force_eq_2[2].set(color=WATER_DENSITY_COLOR)
        force_eq_2[4:].set(color=VOLUME_COLOR)
        force_eq_2.next_to(P_diff_eq_4, UP, aligned_edge=LEFT, buff=0.15)

        self.play(AnimationGroup(
            ReplacementTransform(force_eq[3], force_eq_2[5]),
            FadeOut(force_eq[2], shift=SMALL_BUFF*UP),
            ReplacementTransform(P_diff_eq_4[2:], force_eq_2[2:5]),
            ReplacementTransform(force_eq[:2], force_eq_2[:2]),
            FadeOut(P_diff_eq_4[:2]),
            lag_ratio=0.05, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        indicator = SurroundingRectangle(force_eq_2[-2:], buff=0.05, stroke_width=2)
        sweep = Rectangle(
            width=obj.width, height=0.04,
            stroke_width=0, fill_opacity=1,
            color=VOLUME_COLOR
        ).move_to(obj, aligned_edge=DOWN)

        self.play(Create(indicator))
        self.wait()

        self.play(GrowFromEdge(sweep, LEFT))
        self.play(sweep.animate.stretch_to_fit_height(obj.height, about_edge=DOWN))

        self.wait(1.5)

        force_eq_3 = MathTex(BUOYANT_FORCE_LABEL, '=', WATER_DENSITY_LABEL, 'g', 'V').scale(0.5)
        force_eq_3[0].set(color=BUOYANT_FORCE_COLOR)
        force_eq_3[2].set(color=WATER_DENSITY_COLOR)
        force_eq_3[4].set(color=VOLUME_COLOR)
        force_eq_3.shift(force_eq_2[0].get_center() - force_eq_3[0].get_center())

        self.play(AnimationGroup(
            FadeOut(indicator),
            ReplacementTransform(force_eq_2[4:], force_eq_3[4]),
        ))

        self.remove(*force_eq_2)
        self.add(force_eq_3)

        self.wait(0.2)

        self.play(FadeOut(sweep))

        self.wait(2)

        gravity_vec = force_vec.copy().set(color=GRAVITY_FORCE_COLOR).rotate(PI).next_to(obj, UP, buff=0.05)
        
        gravity_eq = MathTex(GRAVITY_FORCE_LABEL, '=', 'm', 'g')
        gravity_eq.scale(0.5).next_to(force_eq_3, UP, buff=0.15, aligned_edge=LEFT)
        gravity_eq[0].set(color=GRAVITY_FORCE_COLOR)

        shift_val = Group(gravity_eq, force_eq_3).get_y() * DOWN
        gravity_eq.shift(shift_val)

        self.play(AnimationGroup(
            force_eq_3.animate.shift(shift_val).set_opacity(0.5),
            force_vec.animate.set_opacity(0.5),
            FadeIn(gravity_vec, shift=0.25*DOWN),
            FadeIn(gravity_eq, shift=SMALL_BUFF*RIGHT),
            lag_ratio=0.2, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        gravity_eq_2 = MathTex(GRAVITY_FORCE_LABEL, '=', OBJECT_DENSITY_LABEL, 'V', 'g')
        gravity_eq_2.scale(0.5).next_to(force_eq_3, UP, buff=0.15, aligned_edge=LEFT)
        gravity_eq_2[0].set(color=GRAVITY_FORCE_COLOR)
        gravity_eq_2[2].set(color=OBJECT_DENSITY_COLOR)
        gravity_eq_2[3].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            gravity_eq[3].animate.move_to(gravity_eq_2[4]),
            ReplacementTransform(gravity_eq[2], gravity_eq_2[2:4]),
            lag_ratio=0.1, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*gravity_eq)
        self.add(gravity_eq_2)

        gravity_eq_3 = MathTex(GRAVITY_FORCE_LABEL, '=', OBJECT_DENSITY_LABEL, 'g', 'V')
        gravity_eq_3.scale(0.5).next_to(force_eq_3, UP, buff=0.15, aligned_edge=LEFT)
        gravity_eq_3[0].set(color=GRAVITY_FORCE_COLOR)
        gravity_eq_3[2].set(color=OBJECT_DENSITY_COLOR)
        gravity_eq_3[4].set(color=VOLUME_COLOR)

        self.play(AnimationGroup(
            gravity_eq_2[3].animate.move_to(gravity_eq_3[4]),
            gravity_eq_2[4].animate.move_to(gravity_eq_3[3]),
            run_time=1.5, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.remove(*gravity_eq_2)
        self.add(gravity_eq_3)

        self.wait(0.5)

        self.play(AnimationGroup(
            force_vec.animate.set_opacity(1),
            force_eq_3.animate.set_opacity(1),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            self.camera.frame.animate.scale(1.15).shift(0.4*UP),
            Group(gravity_eq_3, force_eq_3).animate.arrange().next_to(glass, UP, buff=0.8),
            lag_ratio=0.1, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        density_equal = MathTex(OBJECT_DENSITY_LABEL, '=', WATER_DENSITY_LABEL).scale(0.5).next_to(glass, UP, buff=SMALL_BUFF)
        density_equal[0].set(color=OBJECT_DENSITY_COLOR)
        density_equal[2].set(color=WATER_DENSITY_COLOR)

        self.play(
            FadeIn(density_equal, shift=SMALL_BUFF*UP),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        left_obj = obj.copy()
        left_force_vec = force_vec.copy()
        left_gravity_vec = gravity_vec.copy()
        left_scene = Group(glass.copy(), water.copy(), left_obj, left_force_vec, left_gravity_vec).shift((glass.width + 1) * LEFT)

        self.play(
            FadeIn(left_scene, target_position=ORIGIN),
            run_time=4
        )

        self.wait()

        density_less = MathTex(OBJECT_DENSITY_LABEL, '<', WATER_DENSITY_LABEL).scale(0.5).next_to(left_scene, UP, buff=SMALL_BUFF)
        density_less[0].set(color=OBJECT_DENSITY_COLOR)
        density_less[2].set(color=WATER_DENSITY_COLOR)

        self.play(
            FadeIn(density_less, shift=SMALL_BUFF*UP),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.5)

        self.play(
            left_gravity_vec.animate.scale(0.67, about_edge=DOWN),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(
            Group(left_obj, left_force_vec, left_gravity_vec).animate.shift((glass.height - obj.height)/2 * UP),
            run_time=5, rate_func=rate_functions.ease_in_quad
        )

        self.wait(2)

        right_obj = obj.copy()
        right_force_vec = force_vec.copy()
        right_gravity_vec = gravity_vec.copy()
        right_scene = Group(glass.copy(), water.copy(), right_obj, right_force_vec, right_gravity_vec).shift((glass.width + 1) * RIGHT)

        self.play(
            FadeIn(right_scene, target_position=ORIGIN),
            run_time=4
        )

        self.wait()

        density_greater = MathTex(OBJECT_DENSITY_LABEL, '>', WATER_DENSITY_LABEL).scale(0.5).next_to(right_scene, UP, buff=SMALL_BUFF)
        density_greater[0].set(color=OBJECT_DENSITY_COLOR)
        density_greater[2].set(color=WATER_DENSITY_COLOR)

        self.play(
            FadeIn(density_greater, shift=SMALL_BUFF*UP),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.5)

        self.play(
            right_gravity_vec.animate.scale(1.5, about_edge=DOWN),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(
            Group(right_obj, right_force_vec, right_gravity_vec).animate.shift(((glass.height - obj.height)/2 - glass_buff) * DOWN),
            run_time=5, rate_func=rate_functions.ease_in_quad
        )

        self.wait(2)
