from manim import *

###########################

POSITION_COLOR = RED
VELOCITY_COLOR = BLUE
ACCELERATION_COLOR = YELLOW

###########################

class CentripetalAccelerationScene(MovingCameraScene):
    def construct(self):
        axes = Axes(x_length=20, y_length=20, tips=False, axis_config={'include_ticks': False, 'color': '#333333'})
        circle = Circle(radius=2, color='#333333', stroke_width=2)
        self.add(axes, circle)

        omega = 0.8
        time = ValueTracker(0)

        dot = Dot(radius=0.1, z_index=100).move_to(2*RIGHT)
        dot.add_updater(lambda mob: mob.move_to(circle.point_at_angle(omega*time.get_value())))

        self.add(dot)

        time.add_updater(lambda tracker, dt: tracker.increment_value(omega*dt))
        self.add(time)

        self.wait(2.5)

        vector_config = {'max_tip_length_to_length_ratio': 0.2, 'max_stroke_width_to_length_ratio': 100}

        v_vector = Vector(color=VELOCITY_COLOR, z_index=50, **vector_config)
        v_vector_size = ValueTracker(0.2)

        def v_vector_updater(mob, dt):
            t = time.get_value()
            dir = -np.sin(omega*t)*RIGHT + np.cos(omega*t)*UP
            mob.put_start_and_end_on(ORIGIN, v_vector_size.get_value() * dir).shift(dot.get_center())
            mob.get_tip().set_opacity(min(1, 4*(v_vector_size.get_value()-0.2)**3))

        v_vector.add_updater(v_vector_updater, call_updater=True)

        self.add(v_vector, v_vector.get_tip())
        self.play(v_vector_size.animate.set_value(1.5), run_time=2, rate_func=rate_functions.ease_out_cubic)

        self.wait(3)

        r_vector = Vector(color=POSITION_COLOR, z_index=50, **vector_config)
        r_vector_size = ValueTracker(0.2)

        def r_vector_updater(mob, dt):
            t = time.get_value()
            dir = np.cos(omega*t)*RIGHT + np.sin(omega*t)*UP
            mob.put_start_and_end_on(ORIGIN, r_vector_size.get_value() * dir)
            mob.get_tip().set_opacity(min(1, 4*(r_vector_size.get_value()-0.2)**3))

        r_vector.add_updater(r_vector_updater, call_updater=True)

        self.add(r_vector)
        self.play(r_vector_size.animate.set_value(2-0.1), run_time=2, rate_func=rate_functions.ease_out_cubic)
        
        self.wait_until(lambda: 0 <= dot.get_y() <= 0.1 and dot.get_x() > 0)

        num_ghosts = 12
        r_ghosts = Group(*[
            Group(
                dot.copy().clear_updaters().move_to(2*RIGHT),
                v_vector.copy().clear_updaters().put_start_and_end_on(2*RIGHT, 2*RIGHT + v_vector_size.get_value()*UP),
                r_vector.copy().clear_updaters().put_start_and_end_on(ORIGIN, r_vector_size.get_value()*RIGHT),
            ).rotate(theta, about_point=ORIGIN)
            for theta in np.linspace(0, 2*PI, num_ghosts, endpoint=False)
        ])

        for ghost in r_ghosts.submobjects:
            for mob in ghost.submobjects:
                mob.set_opacity(0.3)

        for i in range(num_ghosts):
            self.add(r_ghosts[i])
            self.wait(2*PI/(omega**2*num_ghosts))

        time.suspend_updating()

        while time.get_value() >= 2*PI/omega:
            time.increment_value(-2*PI/omega)

        self.wait()

        self.play(AnimationGroup(
            FadeOut(r_vector),
            *[mob.animate.set_opacity(0.08) for ghost in r_ghosts.submobjects for mob in ghost.submobjects],
            self.camera.frame.animate.move_to(dot).scale(0.5),
            time.animate.set_value(0),
            run_time=4.5
        ))

        time.suspend_updating()

        self.wait()

        a_vector = Vector(LEFT, color=ACCELERATION_COLOR, z_index=50, **vector_config).shift(dot.get_center())

        self.play(GrowArrow(a_vector), run_time=2, rate_func=rate_functions.ease_out_cubic)

        self.wait(1.5)

        self.play(
            a_vector.animate.shift(v_vector_size.get_value()*UP),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        def a_vector_updater(mob, dt):
            t = time.get_value()
            dir = -np.cos(omega*t)*RIGHT - np.sin(omega*t)*UP
            mob.put_start_and_end_on(ORIGIN, dir).shift(v_vector.get_end())

        a_vector.add_updater(a_vector_updater)

        v_ghosts = Group(*[
            Group(
                v_vector.copy().clear_updaters().put_start_and_end_on(ORIGIN, v_vector_size.get_value()*UP),
                a_vector.copy().clear_updaters().put_start_and_end_on(ORIGIN, LEFT).shift(v_vector_size.get_value()*UP),
            ).rotate(theta, about_point=ORIGIN) 
            for theta in np.linspace(0, 2*PI, num_ghosts, endpoint=False)
        ])

        for ghost in v_ghosts.submobjects:
            for mob in ghost.submobjects:
                mob.set_opacity(0.3)

        v_ghosts.add_updater(lambda mob: mob.move_to(dot))

        self.camera.frame.add_updater(lambda mob: mob.move_to(dot))

        time.resume_updating()

        for i in range(num_ghosts):
            self.add(v_ghosts[i])
            v_ghosts[i].add_updater(lambda mob: mob.shift(dot.get_center() - mob[0].get_start()))
            self.wait(2*PI/(omega**2*num_ghosts))

        time.clear_updaters()
        self.camera.frame.clear_updaters()

        time.set_value(0)

        self.wait(2)

        dot.clear_updaters()
        r_vector.clear_updaters()
        v_vector.clear_updaters()
        a_vector.clear_updaters()

        dot_copy = dot.copy()
        v_vector_copy = v_vector.copy()
        self.add(dot_copy, v_vector_copy)

        for ghost in v_ghosts:
            ghost.clear_updaters()

        self.play(AnimationGroup(
            FadeIn(r_vector.put_start_and_end_on(ORIGIN, r_vector_size.get_value()*RIGHT)),
            *[mob.animate.set_opacity(0.3) for ghost in r_ghosts.submobjects for mob in ghost.submobjects],
            self.camera.frame.animate.scale(2).shift(RIGHT),
            FadeOut(axes),
            Group(dot_copy, v_vector_copy, a_vector, *v_ghosts).animate.shift(4*RIGHT),
            run_time=4.5
        ))

        self.wait()

        r_circle = circle.copy().set(stroke_width=4, color=WHITE)
        v_circle = Circle(radius=v_vector_size.get_value(), color=WHITE).move_to(v_vector_copy.get_start()).rotate(90*DEGREES)

        self.play(AnimationGroup(
            Rotate(Group(dot, r_vector, v_vector), 2*PI, about_point=r_circle.get_center()),
            Create(r_circle),
            Rotate(Group(v_vector_copy, a_vector), 2*PI, about_point=v_circle.get_center()),
            Create(v_circle),
            run_time=2*PI/(omega**2), rate_func=rate_functions.ease_in_out_sine
        ))

        self.wait(2)

        r_T_eq = MathTex('T', '=', '{2', r'\pi', 'r', r'\over', 'v}').scale(0.9).next_to(r_circle, DOWN, buff=0.8)
        r_T_eq[4].set(color=POSITION_COLOR)
        r_T_eq[6].set(color=VELOCITY_COLOR)

        self.play(AnimationGroup(
            self.camera.frame.animate.set_y(Group(r_T_eq, r_circle).get_y()),
            FadeIn(Group(r_T_eq[:2], r_T_eq[5]), shift=0.5*DOWN),
            lag_ratio=0.05, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            ReplacementTransform(r_circle, r_T_eq[2:4].copy().set_fill(opacity=0)),
            FadeIn(r_T_eq[2:4]),
            ReplacementTransform(r_vector, r_T_eq[4]),
            lag_ratio=0.3, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.play(
            ReplacementTransform(v_vector, r_T_eq[6]),
            run_time=2, rate_func=rate_functions.ease_in_out_cubic
        )

        self.play(FadeOut(dot))

        self.wait(2)

        v_T_eq = MathTex('T', '=', '{2', r'\pi', 'v', r'\over', 'a}').scale(0.9).set_x(v_circle.get_x()).set_y(r_T_eq.get_y())
        v_T_eq[4].set(color=VELOCITY_COLOR)
        v_T_eq[6].set(color=ACCELERATION_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(r_T_eq[:2], v_T_eq[:2]),
            FadeIn(v_T_eq[5], shift=0.2*RIGHT),
            lag_ratio=0.3, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.play(AnimationGroup(
            ReplacementTransform(v_circle, v_T_eq[2:4].copy().set_fill(opacity=0)),
            FadeIn(v_T_eq[2:4]),
            ReplacementTransform(v_vector_copy, v_T_eq[4]),
            lag_ratio=0.3, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.play(
            ReplacementTransform(a_vector, v_T_eq[6]),
            run_time=2, rate_func=rate_functions.ease_in_out_cubic
        )

        self.play(dot_copy.animate.set_fill(opacity=0.3))

        self.wait(3)

        eq_1 = MathTex('{2', r'\pi', 'r', r'\over', 'v}', '=', '{2', r'\pi', 'v', r'\over', 'a}')
        eq_1.set_x(3).set_y(r_T_eq.get_y())
        eq_1[2].set(color=POSITION_COLOR)
        eq_1[4:9:4].set(color=VELOCITY_COLOR)
        eq_1[10].set(color=ACCELERATION_COLOR)

        self.play(AnimationGroup(
            FadeOut(v_T_eq[0]),
            ReplacementTransform(r_T_eq[2:], eq_1[:5]),
            ReplacementTransform(v_T_eq[1:], eq_1[5:]),
            self.camera.frame.animate.move_to(eq_1).scale(0.8),
            lag_ratio=0.2, run_time=5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(
            FadeOut(Group(eq_1[:2], eq_1[6:8]), shift=0.5*UP),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.5)

        eq_2 = MathTex('a', '=', '{v', '^2', r'\over', 'r}').move_to(eq_1)
        eq_2[0].set(color=ACCELERATION_COLOR)
        eq_2[2].set(color=VELOCITY_COLOR)
        eq_2[5].set(color=POSITION_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(eq_1[10], eq_2[0]),
            ReplacementTransform(eq_1[8], eq_2[2]),
            ReplacementTransform(eq_1[4], eq_2[2]),
            ReplacementTransform(eq_1[2], eq_2[5]),
            ReplacementTransform(eq_1[5], eq_2[1]),
            ReplacementTransform(eq_1[3], eq_2[4]),
            ReplacementTransform(eq_1[9], eq_2[4]),
            FadeIn(eq_2[3], target_position=eq_2[2]),
            lag_ratio=0.1, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()
