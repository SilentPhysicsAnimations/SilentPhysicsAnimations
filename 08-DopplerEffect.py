from manim import *

########################

FREQUENCY_COLOR = YELLOW
WAVELENGTH_COLOR = BLUE
VELOCITY_COLOR = RED

SOUND_SPEED_LABEL = 'c'
SOUND_COLOR = WHITE

########################

TEX = { 'tex_template': TexTemplateLibrary.default }

########################

class Flasher(Mobject):
    def __init__(self, mobject, flash_color=YELLOW, flash_speed=1, **kwargs):
        super().__init__(**kwargs)
        self.add(mobject)
        self.mobject = mobject
        self.base_color = mobject.color
        self.flash_color = flash_color
        self.flash_speed = flash_speed
        self.flash_timer = ValueTracker(0)
        self._flash_updater = None

    def proceed_flash(self, dt):
        alpha = self.flash_timer.get_value()
        self.mobject.set(color=color.interpolate_color(
            self.base_color,
            self.flash_color,
            alpha
        ))
        self.flash_timer.set_value(max(0, alpha - self.flash_speed*dt))
        if self.flash_timer.get_value() == 0:
            self.remove_updater(self._flash_updater)
            self._flash_updater = None

    def flash(self):
        self.flash_timer.set_value(1)
        if self._flash_updater is None:
            self._flash_updater = lambda mob, dt: mob.proceed_flash(dt)
            self.add_updater(self._flash_updater)


class DopplerEffectScene(MovingCameraScene):
    def construct(self):
        c = 1.33
        f = 1.5

        emitter = Flasher(Dot(radius=0.5, z_index=100), flash_speed=c*f, flash_color=FREQUENCY_COLOR)
        self.add(emitter)
        lwave = Arc(
            radius=0.5*np.sqrt(2)-0.02,
            start_angle=135*DEGREES, angle=90*DEGREES,
            stroke_color=SOUND_COLOR, fill_color=BLACK
        ).move_to(ORIGIN, aligned_edge=RIGHT)
        rwave = Arc(
            radius=0.5*np.sqrt(2)-0.02,
            start_angle=-45*DEGREES, angle=90*DEGREES,
            stroke_color=SOUND_COLOR, fill_color=BLACK
        ).move_to(ORIGIN, aligned_edge=LEFT)

        listener = Flasher(Dot(radius=0.5, z_index=100), flash_speed=c*f, flash_color=FREQUENCY_COLOR)
        listener.shift(12*RIGHT)
        self.add(listener)

        listener_should_flash = False
        def wave_updater(dir):
            def updater(mob, dt):
                mob.shift(c*dt*dir)
                mob.set_opacity((1-0.1*dt)*mob.get_stroke_opacity())
                if mob.get_x() >= listener.get_x() and listener_should_flash:
                    listener.flash()
                if mob.get_stroke_opacity() <= 0.01 or mob.get_x() >= listener.get_x():
                    mob.set_opacity(0)
                    mob.clear_updaters()
                    self.remove(mob)
            return updater
        
        wave_time = 0
        emitter_should_flash = False
        def scene_updater(dt):
            nonlocal wave_time
            wave_time -= dt
            if wave_time <= 0.85/f and wave_time + dt > 0.85/f and emitter_should_flash:
                emitter.flash()
            if wave_time <= 0:
                wave_time += 1/f
                self.add(
                    lwave.copy().shift(emitter.get_center()).add_updater(wave_updater(LEFT)),
                    rwave.copy().shift(emitter.get_center()).add_updater(wave_updater(RIGHT))
                )

        self.add_updater(scene_updater)

        self.wait(5)
        self.play(AnimationGroup(
            self.camera.frame.animate.shift(4*RIGHT),
            listener.animate.move_to(8*RIGHT),
            run_time=3, lag_ratio=0.4
        ))
        self.wait(1)
        listener_should_flash = True
        self.wait(2)

        f_prime_label = MathTex('f', "'", color=FREQUENCY_COLOR).scale(0.9).next_to(listener, DOWN)
        f_prime_label.shift(f_prime_label[0].get_top() - f_prime_label.get_top())
        
        self.play(
            FadeIn(f_prime_label[0], target_position=listener),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(1.5)
        listener_should_flash = False
        self.wait(1)
        emitter_should_flash = True
        self.wait(2)

        f_label = MathTex('f', color=FREQUENCY_COLOR).scale(0.9).move_to(f_prime_label[0])
        self.add(f_label)
        self.play(
            f_label.animate.next_to(emitter, DOWN),
            run_time=2, rate_func=rate_functions.ease_in_out_cubic
        )

        self.wait(1.5)
        emitter_should_flash = False
        self.wait(2)
        self.wait_until(lambda: wave_time <= 1/config.frame_rate + 1e-3)

        wavelength = c/f
        lambda_line = Line(ORIGIN, wavelength*RIGHT, color=WAVELENGTH_COLOR, z_index=50)
        lambda_line.shift(rwave.width/2*RIGHT)
        self.play(
            GrowFromEdge(lambda_line, LEFT),
            run_time=1/f, rate_func=rate_functions.linear
        )

        lambda_line.add_updater(lambda mob, dt: mob.shift(c*dt*RIGHT), call_updater=True)
        self.wait_until(lambda: lambda_line.get_right()[0] >= 3.5)
        lambda_line.clear_updaters()
        
        self.play(
            lambda_line.animate.move_to(4*RIGHT + 2*DOWN),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )
        
        self.wait(0.5)
        
        lambda_label = MathTex(r'\lambda', color=WAVELENGTH_COLOR).scale(0.9).next_to(lambda_line, UP)

        self.play(
            FadeIn(lambda_label, target_position=lambda_line),
            run_time=1.5, rate_funct=rate_functions.ease_out_cubic
        )

        self.wait(1)

        v_vector = Vector(RIGHT, color=VELOCITY_COLOR).next_to(emitter, UP).shift(0.5*RIGHT)
        v_label = MathTex('v', color=VELOCITY_COLOR).scale(0.9).next_to(v_vector, UP)

        self.play(AnimationGroup(
            FadeIn(v_vector, shift=UP),
            FadeIn(v_label, shift=UP),
            lag_ratio=0.2, run_time=1.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v = c/3
        for mob in [emitter, f_label, v_vector, v_label]:
            mob.add_updater(lambda mob, dt: mob.shift(v*dt*RIGHT))

        self.wait()
        self.wait_until(lambda: wave_time <= 1/config.frame_rate + 1e-3)

        wavelength_prime = (c - v)/f

        lambda_prime_line = Line(ORIGIN, wavelength_prime*RIGHT, color=WAVELENGTH_COLOR, z_index=50)
        lambda_prime_line.shift(rwave.width*RIGHT)
        self.add(lambda_prime_line)

        lambda_prime_line.add_updater(lambda mob, dt: mob.shift(c*dt*RIGHT), call_updater=True)
        self.wait_until(lambda: lambda_prime_line.get_right()[0] >= 4)
        lambda_prime_line.clear_updaters()
        
        self.play(
            lambda_prime_line.animate.next_to(lambda_line, DOWN, aligned_edge=LEFT),
            run_time=1.4, rate_func=rate_functions.ease_out_cubic
        )
        
        self.wait(0.3)
        
        lambda_prime_label = MathTex(r"\lambda'", color=WAVELENGTH_COLOR).scale(0.9).next_to(lambda_prime_line, DOWN).set_x(4)

        self.play(
            FadeIn(lambda_prime_label, shift=0.25*DOWN),
            run_time=0.8, rate_funct=rate_functions.ease_out_cubic
        )

        self.wait(1.2)

        listener_should_flash = True
        self.wait(1)
        self.play(
            FadeIn(f_prime_label[1], target_position=f_prime_label[0]),
            rate_func=rate_functions.ease_out_cubic
        )
        self.wait(2)

        for mob in [v_vector, v_label, f_label]:
            mob.clear_updaters()

        scene_shift = 2*DOWN
        v_vector_shift_val = lambda_prime_line.get_right() - v_vector.get_left() + (v/f + 0.1)*RIGHT

        self.play(AnimationGroup(
            self.camera.frame.animate.move_to(Group(lambda_line, lambda_prime_line)).shift(scene_shift).scale(0.4),
            Group(lambda_line, lambda_prime_line).animate.shift(scene_shift),
            lambda_label.animate.scale(0.4).next_to(lambda_line, UP, buff=0.1).shift(scene_shift),
            lambda_prime_label.animate.scale(0.4).next_to(lambda_prime_line, DOWN, buff=0.1).shift(scene_shift),
            f_label.animate.scale(0.4).next_to(lambda_line, LEFT, buff=0.1).shift(scene_shift),
            f_prime_label.animate.scale(0.4).next_to(lambda_prime_line, LEFT, buff=0.1).shift(scene_shift),
            v_vector.animate.shift(v_vector_shift_val + scene_shift).set(stroke_width=4),
            v_label.animate.shift(v_vector_shift_val + scene_shift + 0.25*DOWN).scale(0.4),
            run_time=3.5
        ))

        emitter.clear_updaters()
        listener_should_flash = False

        self.wait(2)

        lambda_diff_label = MathTex('v', '/', 'f').scale(0.4*0.9).move_to(v_label)
        lambda_diff_label[0].set(color=VELOCITY_COLOR)
        lambda_diff_label[2].set(color=FREQUENCY_COLOR)

        v_line = Line(v_vector.get_start(), v_vector.get_end(), color=VELOCITY_COLOR)
        self.add(v_line)

        self.play(AnimationGroup(
            FadeOut(v_vector, rate_func=rate_functions.ease_out_cubic),
            v_label.animate(rate_func=rate_functions.ease_out_cubic).move_to(lambda_diff_label[0]),
            FadeIn(lambda_diff_label[1:], shift=0.1*LEFT, rate_func=rate_functions.ease_out_cubic),
            v_line.animate(rate_func=rate_functions.ease_out_cubic).stretch_to_fit_width(v/f),
            lag_ratio=0.4, run_time=3
        ))

        self.remove(v_label)
        self.add(lambda_diff_label)

        self.wait(1.5)
        
        self.play(
            v_line.animate.next_to(lambda_prime_line, buff=0),
            lambda_diff_label.animate.move_to(lambda_prime_label).shift(wavelength/2*RIGHT),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_in_out_cubic,
        )

        self.wait(1)

        lambda_eq = MathTex(r'\lambda', '=', r"\lambda'", '+', '{v', r'\over', 'f}').scale(0.4*0.9)
        lambda_eq[0:3:2].set(color=WAVELENGTH_COLOR)
        lambda_eq[4].set(color=VELOCITY_COLOR)
        lambda_eq[6].set(color=FREQUENCY_COLOR)
        lambda_eq.next_to(Group(lambda_line, lambda_prime_line))

        self.play(AnimationGroup(
            ReplacementTransform(lambda_label.copy(), lambda_eq[0]),
            FadeIn(lambda_eq[1]),
            ReplacementTransform(lambda_prime_label.copy(), lambda_eq[2]),
            FadeIn(lambda_eq[3]),
            ReplacementTransform(lambda_diff_label[0], lambda_eq[4]),
            ReplacementTransform(lambda_diff_label[1], lambda_eq[5]),
            ReplacementTransform(lambda_diff_label[2], lambda_eq[6]),
            FadeOut(v_line),
            lag_ratio=0.05, run_time=4,
        ))

        self.wait(1.5)

        c = 0.4
        rwave_c = rwave.copy().scale(0.4).set(stroke_width=2).move_to(
            mid(self.camera.frame.get_left(), self.camera.frame.get_corner(UL)),
            aligned_edge=RIGHT
        ).shift(0.5*LEFT)
        c_vector = Vector(0.33*RIGHT, max_tip_length_to_length_ratio=0.2, color=SOUND_COLOR).next_to(rwave_c, buff=0)
        rwave_c.add_updater(wave_updater(RIGHT))
        c_vector.add_updater(wave_updater(RIGHT))
        c_label = MathTex(SOUND_SPEED_LABEL, color=SOUND_COLOR, **TEX).scale(0.4*0.9)
        c_label.add_updater(lambda mob: mob.next_to(c_vector, UP, buff=0.1), call_updater=True)

        self.add(rwave_c, c_vector, c_label)

        self.wait_until(lambda: rwave_c.get_x() >= 3.5)
        c_label.clear_updaters()

        lambda_formula = MathTex(r'\lambda', '=', '{' + SOUND_SPEED_LABEL, r'\over', 'f}', **TEX)
        lambda_formula.scale(0.4*0.9).move_to(lambda_label, aligned_edge=DOWN)

        lambda_prime_formula = MathTex(r"\lambda'", '=', '{' + SOUND_SPEED_LABEL, r'\over', "f'}", **TEX)
        lambda_prime_formula.scale(0.4*0.9).move_to(lambda_prime_label, aligned_edge=UP)

        lambda_formula[0].set(color=WAVELENGTH_COLOR)
        lambda_formula[2].set(color=SOUND_COLOR)
        lambda_formula[4].set(color=FREQUENCY_COLOR)
        lambda_prime_formula[0].set(color=WAVELENGTH_COLOR)
        lambda_prime_formula[2].set(color=SOUND_COLOR)
        lambda_prime_formula[4].set(color=FREQUENCY_COLOR)

        c_label_copy = c_label.copy()

        self.play(AnimationGroup(
            c_label_copy.animate.move_to(lambda_formula[2]),
            c_label.animate.move_to(lambda_prime_formula[2]),

            lambda_label.animate.move_to(lambda_formula[0]),
            FadeIn(lambda_formula[1]),
            FadeIn(lambda_formula[3]),
            f_label.animate.move_to(lambda_formula[4]),

            lambda_prime_label.animate.move_to(lambda_prime_formula[0]),
            FadeIn(lambda_prime_formula[1]),
            FadeIn(lambda_prime_formula[3]),
            f_prime_label.animate.move_to(lambda_prime_formula[4]),

            lag_ratio=0.05, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(c_label, c_label_copy, f_label, f_prime_label, lambda_label, lambda_prime_label)
        self.add(lambda_formula, lambda_prime_formula)

        self.wait()

        no_lambda_eq = MathTex(
            '{' + SOUND_SPEED_LABEL, r'\over', 'f}', '=', '{' + SOUND_SPEED_LABEL, '\over', "f'}", '+', '{v', r'\over', 'f}',
            **TEX
        ).scale(0.4*0.9)
        no_lambda_eq[0:5:4].set(color=SOUND_COLOR)
        no_lambda_eq[2::4].set(color=FREQUENCY_COLOR)
        no_lambda_eq[8].set(color=VELOCITY_COLOR)
        no_lambda_eq.move_to(lambda_eq)

        self.play(AnimationGroup(
            FadeOut(lambda_line),
            FadeOut(lambda_prime_line),

            FadeOut(lambda_eq[0], shift=0.1*DOWN),
            FadeOut(lambda_eq[2], shift=0.1*DOWN),

            FadeOut(lambda_formula[0:2]),
            AnimationGroup(
                lambda_formula[2].animate.move_to(no_lambda_eq[0]),
                lambda_formula[3].animate.move_to(no_lambda_eq[1]),
                lambda_formula[4].animate.move_to(no_lambda_eq[2]),
            ),

            FadeOut(lambda_prime_formula[0:2]),
            AnimationGroup(
                lambda_prime_formula[2].animate.move_to(no_lambda_eq[4]),
                lambda_prime_formula[3].animate.move_to(no_lambda_eq[5]),
                lambda_prime_formula[4].animate.move_to(no_lambda_eq[6]),
            ),

            AnimationGroup(
                lambda_eq[1].animate.move_to(no_lambda_eq[3]),
                lambda_eq[3].animate.move_to(no_lambda_eq[7]),
                lambda_eq[4].animate.move_to(no_lambda_eq[8]),
                lambda_eq[5].animate.move_to(no_lambda_eq[9]),
                lambda_eq[6].animate.move_to(no_lambda_eq[10]),
            ),

            self.camera.frame.animate.move_to(no_lambda_eq).scale(0.8),

            lag_ratio=0.05, run_time=4
        ))

        self.remove(*lambda_formula, *lambda_prime_formula, *lambda_eq)
        self.add(no_lambda_eq)
        
        self.wait(2)

        eq_2 = MathTex('{' + SOUND_SPEED_LABEL, '-', 'v', r'\over', 'f}', '=', '{' + SOUND_SPEED_LABEL, r'\over', "f'}", **TEX)
        eq_2.scale(0.4*0.9).move_to(no_lambda_eq)
        eq_2[0::6].set(color=SOUND_COLOR)
        eq_2[4::4].set(color=FREQUENCY_COLOR)
        eq_2[2].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            *[no_lambda_eq[i].animate.move_to(eq_2[j])
                for i, j in [(0, 0), (2, 4), (3, 5), (4, 6), (5, 7), (6, 8), (8, 2), (10, 4)]],
            ReplacementTransform(no_lambda_eq[1], eq_2[3]),
            ReplacementTransform(no_lambda_eq[7], eq_2[1]),
            FadeOut(no_lambda_eq[9]),
            run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*no_lambda_eq)
        self.add(eq_2)

        self.wait(1.5)

        eq_3 = MathTex('{1', '-', '{v', r'\over', SOUND_SPEED_LABEL + '}', r'\over', 'f}', '=', '{1', r'\over', "f'}", **TEX)
        eq_3.scale(0.4*0.9)
        eq_3[4].set(color=SOUND_COLOR)
        eq_3[6::4].set(color=FREQUENCY_COLOR)
        eq_3[2].set(color=VELOCITY_COLOR)
        eq_3.shift(eq_2[5].get_center() - eq_3[7].get_center())

        self.play(AnimationGroup(
            *[eq_2[i].animate.move_to(eq_3[j])
                for i, j in [(1, 1), (4, 6), (5, 7), (6, 4), (7, 9), (8, 10)]],
            ReplacementTransform(eq_2[0], eq_3[0]),
            ReplacementTransform(eq_2[2], eq_3[2]),
            ReplacementTransform(eq_2[3], eq_3[5]),
            ReplacementTransform(eq_2[6], eq_3[4]),
            FadeIn(eq_3[3]),
            FadeIn(eq_3[8]),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(*eq_2)
        self.add(eq_3)

        self.wait(1.5)

        final_eq = MathTex("f'", '=', '{1', r'\over', '1', '-', '{v', r'\over', SOUND_SPEED_LABEL + '}}', 'f', **TEX)
        final_eq.scale(0.4*0.9)
        final_eq[8].set(color=SOUND_COLOR)
        final_eq[0::9].set(color=FREQUENCY_COLOR)
        final_eq[6].set(color=VELOCITY_COLOR)
        final_eq.move_to(self.camera.frame)

        self.play(AnimationGroup(
            *[eq_3[i].animate.move_to(final_eq[j])
                for i, j in [(0, 4), (1, 5), (2, 6), (3, 7), (4, 8), (6, 9), (7, 1), (8, 2)]],
            FadeOut(eq_3[5]),
            ReplacementTransform(eq_3[9], final_eq[3]),
            ReplacementTransform(eq_3[10], final_eq[0]),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))
        
        self.remove(*eq_3)
        self.add(final_eq)

        self.wait(2)
