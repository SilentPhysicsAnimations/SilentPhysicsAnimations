from manim import *

####################
 
LIGHT_COLOR = YELLOW
TIME_COLOR = BLUE
VELOCITY_COLOR = RED

####################

class TrainCar(VMobject):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.interior = Rectangle(height=1, width=16/9, fill_color=BLACK, fill_opacity=1, **kwargs)
        self.left_wheel = Dot(radius=0.25, fill_color=GRAY).move_to(self.interior.get_corner(DL), aligned_edge=LEFT)
        self.right_wheel = Dot(radius=0.25, fill_color=GRAY).move_to(self.interior.get_corner(DR), aligned_edge=RIGHT)
        self.add(self.left_wheel, self.right_wheel, self.interior)


class TimeDilationScene(MovingCameraScene):
    def construct(self):
        car = TrainCar().scale(1.5)

        car_buff = 0.6
        other_cars = Group(*[car.copy() for _ in range(3+7)])
        train = Group(*other_cars[:3], car, *other_cars[3:]).arrange(buff=car_buff)
        connector = Line(ORIGIN, (car_buff-0.04)*RIGHT, stroke_width=8)
        connectors = Group(*[
            connector.copy().move_to(mid(left_car.get_right(), left_car.get_corner(DR)), aligned_edge=LEFT).shift(0.02*RIGHT)
            for left_car in [car, *other_cars[:-1]]
        ])
        train.add(*connectors)

        ground = Rectangle(
            width=train.width, height=0.02,
            stroke_width=0, fill_opacity=1,
        ).move_to(train.get_bottom(), aligned_edge=UP)

        tree_num = 33
        tree_buff = ground.width / (tree_num - 1)
        tree = Line(
            ORIGIN, (self.camera.frame.get_top()[1] - ground.get_top()[1])*UP,
            stroke_width=5, stroke_opacity=0.2,
            sheen_factor=-0.5, sheen_direction=UP
        )
        trees = Group(*[tree.copy() for _ in range(tree_num)]).arrange(buff=tree_buff)
        trees.move_to(ground.get_top(), aligned_edge=DOWN)
        self.add(trees, train, ground)

        dist = car.width + car_buff
        v = 1

        self.camera.frame.move_to(car)

        for _ in range(2):
            self.play(
                train.animate.shift(dist*RIGHT),
                run_time=dist/v, rate_func=rate_functions.linear
            )
            train.shift(dist*LEFT)

        def trees_updater(mob, dt):
            mob.shift(v*dt*LEFT)
            if mob.get_x() - ground.get_x() <= -tree_buff:
                mob.shift(tree_buff*RIGHT)

        trees.add_updater(trees_updater)

        light = Rectangle(
            height=0.04, width=car.interior.width/4,
            stroke_width=0, fill_color=LIGHT_COLOR, fill_opacity=1
        )
        light.move_to(car.interior, aligned_edge=UP).shift(light.height/2 * DOWN)

        glow = Rectangle(
            height=car.interior.height-0.04, width=car.interior.width-0.04,
            stroke_width=0, fill_color=LIGHT_COLOR, fill_opacity=0.12,
            sheen_factor=-1, sheen_direction=DOWN
        )
        glow.move_to(light, aligned_edge=UP)

        self.play(
            AnimationGroup(
                self.camera.frame.animate.move_to(car.interior).scale(0.5),
                *[wheel.animate.set(color='#222222') for wheel in [car.left_wheel for car in other_cars]],
                *[wheel.animate.set(color='#222222') for wheel in [car.right_wheel for car in other_cars]],
                *[interior.animate.set(stroke_color='#444444') for interior in [car.interior for car in other_cars]],
                *[connector.animate.set(color='#444444') for connector in connectors],
                FadeIn(light, target_position=car.get_top()),
                FadeIn(glow)
            ),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        rays = Group(*[Line(ORIGIN, 2*light.height*DOWN, color=LIGHT_COLOR, stroke_width=0.8) for _ in range(6)])
        rays.arrange(buff=light.width/len(rays)).next_to(light, DOWN, buff=light.height)

        self.play(
            *[ShowPassingFlash(ray, time_width=0.5) for ray in rays],
            run_time=0.5, rate_func=rate_functions.linear
        )

        self.wait(1.5)

        stopwatch_frame = Circle(radius=0.2, color=WHITE, fill_opacity=1, fill_color=BLACK, stroke_width=3)
        stopwatch_mark = Line(ORIGIN, 0.185*UP, stroke_width=2)
        stopwatch_hand = Line(ORIGIN, 0.185*UP, stroke_width=2)
        stopwatch_handle = Line(ORIGIN, 0.2*RIGHT, stroke_width=3).next_to(stopwatch_frame, UP, buff=0.065)
        stopwatch_button = Line(0.225*UP, 0.26*UP, stroke_width=3).rotate(-45*DEGREES, about_point=ORIGIN)
        stopwatch = Group(stopwatch_frame, stopwatch_mark, stopwatch_hand, stopwatch_handle, stopwatch_button)
        stopwatch.scale(0.8).move_to(mid(light.get_corner(UR), glow.get_corner(DR))).shift(glow.height/6 * UP)


        self.play(
            FadeIn(stopwatch, shift=0.1*UP),
            run_time=0.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        h = car.interior.height - 2*light.height
        c = 1.5

        ray = Dot(radius=light.height/2, color=LIGHT_COLOR).move_to(light)
        ray_trace = TracedPath(
            ray.get_center,
            stroke_width=4, stroke_color=LIGHT_COLOR, stroke_opacity=[1, 0],
            dissipating_time=(h/c)/6
        )

        self.add(ray, ray_trace)

        stopwatch_sweep = AnnularSector()
        stopwatch_sweep.add_updater(lambda mob: mob.become(
            AnnularSector(
                inner_radius=0, outer_radius=0.145,
                start_angle=90*DEGREES,
                angle=-angle_between_vectors(UP, stopwatch_hand.get_end() - stopwatch_hand.get_start()),
                fill_opacity=0.5, color=TIME_COLOR
            ).shift(stopwatch_frame.get_center())
        ))

        self.add(stopwatch_sweep)

        self.play(
            ray.animate(run_time=h/c, rate_func=rate_functions.linear).shift(h*DOWN),
            Rotate(
                stopwatch_hand, -120*DEGREES,
                about_point=stopwatch_frame.get_center(),
                run_time=h/c, rate_func=rate_functions.linear
            )
        )

        stopwatch_sweep.clear_updaters()

        self.play(
            ray.animate.stretch_to_fit_height(1e-6, about_point=ray.get_bottom()),
            run_time=ray.height/c, rate_func=rate_functions.linear
        )

        self.wait(2)

        self.remove(ray, ray_trace)

        time_A_label = MathTex('T_A', color=TIME_COLOR).scale(0.4)
        time_A_label.move_to(mid(light.get_corner(UR), glow.get_corner(DR))).shift(glow.height/6 * DOWN)
        self.play(ReplacementTransform(stopwatch_sweep, time_A_label))

        self.wait()

        brace_height = glow.height - light.height
        h_brace = BraceBetweenPoints(ORIGIN, 2*brace_height*DOWN, buff=2*light.height)
        h_brace.scale(0.5, about_point=ORIGIN).shift(light.get_corner(DL))
        h_brace_label = MathTex('h').scale(0.4).next_to(h_brace, LEFT, buff=SMALL_BUFF)

        c_vector = Vector(glow.height/3 * DOWN, color=LIGHT_COLOR).next_to(light, DOWN, buff=light.height)
        c_vector_label = MathTex('c', color=LIGHT_COLOR).scale(0.4).next_to(c_vector, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            AnimationGroup(
                GrowArrow(c_vector),
                FadeIn(c_vector_label, shift=SMALL_BUFF*DOWN),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromEdge(h_brace, UP),
                FadeIn(h_brace_label, target_position=h_brace.get_left()),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.4, run_time=4
        ))

        self.wait(2)

        time_A_eq = MathTex('T_A', '=', '{h', r'\over', 'c}').scale(0.4).move_to(time_A_label)
        time_A_eq[0].set(color=TIME_COLOR)
        time_A_eq[4].set(color=LIGHT_COLOR)

        self.play(AnimationGroup(
            AnimationGroup(
                time_A_label.animate.move_to(time_A_eq[0]),
                FadeIn(time_A_eq[1], shift=0.05*RIGHT),
                h_brace_label.animate.move_to(time_A_eq[2]),
                FadeIn(time_A_eq[3]),
                c_vector_label.animate.move_to(time_A_eq[4]),
                lag_ratio=0.05, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeOut(h_brace), FadeOut(c_vector),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=4
        ))

        self.remove(time_A_label, h_brace_label, c_vector_label)
        self.add(*time_A_eq)

        self.wait()

        trees.clear_updaters()

        trees.shift(tree_buff*RIGHT)
        
        tree_shift_dist = trees.get_center()[0] - ground.get_center()[0]
        self.play(
            trees.animate.shift(tree_shift_dist*LEFT),
            run_time=tree_shift_dist/v, rate_func=rate_functions.linear
        )

        train.add(light, glow)
        self.add(train)
        train.add_updater(lambda mob, dt: mob.shift(v*dt*RIGHT))

        self.play(AnimationGroup(
            time_A_eq.animate.scale(1.2).set_x(-0.8*(7-1/9)+0.5, LEFT).set_y(mid(ground.get_y(), -0.8*4)),
            stopwatch.animate.scale(1.8).set_x(0).set_y(mid(ground.get_y(), -0.8*4)),
            self.camera.frame.animate.scale(1.6).center(),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.play(
            Rotate(
                stopwatch_hand, 120*DEGREES,
                about_point=stopwatch_frame.get_center()
            ),
            run_time=0.5
        )

        c = v * np.sqrt(h**2/dist**2 + 1)
        gamma = 1 / np.sqrt(1 - v**2/c**2)
        T_B = gamma * h/c

        self.wait_until(lambda: car.get_x() + v/60 >= -v*T_B/2)
        train.shift((-v*T_B/2 - car.get_x())*LEFT)

        def trace_updater(mob):
            mob.set_points([np.array([ray.get_x(), p[1], 0]) for p in mob.points])

        ray = Dot(radius=light.height/2, color=LIGHT_COLOR).move_to(light)
        ray_trace = TracedPath(
            ray.get_center,
            stroke_width=4, stroke_color=LIGHT_COLOR, stroke_opacity=[1, 0],
            dissipating_time=gamma*(h/c)/6
        ).add_updater(trace_updater)

        path_trace = TracedPath(ray.get_center, stroke_width=3, stroke_opacity=0.4)

        self.add(ray, ray_trace, path_trace)

        stopwatch_sweep = AnnularSector()
        stopwatch_sweep.add_updater(lambda mob: mob.become(
            AnnularSector(
                inner_radius=0, outer_radius=1.8*0.15,
                start_angle=90*DEGREES,
                angle=-angle_between_vectors(UP, stopwatch_hand.get_end() - stopwatch_hand.get_start())
                    if stopwatch_hand.get_x() >= stopwatch_frame.get_x()
                    else -2*PI+angle_between_vectors(UP, stopwatch_hand.get_end() - stopwatch_hand.get_start()),
                fill_opacity=0.5, color=TIME_COLOR
            ).shift(stopwatch_frame.get_center())
        ))

        self.add(stopwatch_sweep)

        self.play(
            ray.animate(run_time=gamma*h/c, rate_func=rate_functions.linear).shift(h*DOWN + v*T_B*RIGHT),
            Rotate(
                stopwatch_hand, -220*DEGREES,
                about_point=stopwatch_frame.get_center(),
                run_time=T_B, rate_func=rate_functions.linear
            )
        )

        path_trace.clear_updaters()
        self.remove(path_trace)
        path_trace_line = Line(v*T_B*LEFT, h*DOWN, stroke_width=3, stroke_opacity=0.4, z_index=1).shift(light.get_center())
        self.add(path_trace_line)

        stopwatch_sweep.clear_updaters()
        train.clear_updaters()

        ray_trace.clear_updaters()

        self.play(
            ray.animate.stretch_to_fit_height(1e-6, about_point=ray.get_bottom()),
            run_time=gamma*ray.height/c, rate_func=rate_functions.linear
        )

        self.remove(ray)

        self.wait()

        time_B_label = MathTex('T_B', color=TIME_COLOR).scale(0.64).next_to(stopwatch)
        self.play(ReplacementTransform(stopwatch_sweep, time_B_label))

        self.wait(2)

        dist_brace = BraceBetweenPoints(ORIGIN, 1.25*v*T_B*RIGHT).scale(0.8)
        dist_brace.next_to(ground, UP, buff=SMALL_BUFF)
        dist_brace.stretch_to_fit_width(1e-6, about_edge=RIGHT)

        self.play(AnimationGroup(
            train.animate.shift(v*T_B*LEFT),
            dist_brace.animate.stretch_to_fit_width(v*T_B, about_edge=RIGHT),
            path_trace_line.animate.set_opacity(0.67).set(color=LIGHT_COLOR),
            Rotate(stopwatch_hand, 220*DEGREES, about_point=stopwatch_frame.get_center()),
            run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        v_vector = Vector(color=VELOCITY_COLOR).shift(mid(car.interior.get_right(), car.interior.get_corner(UR)))
        v_vector_label = MathTex('v', color=VELOCITY_COLOR).scale(0.64).next_to(v_vector.get_tip(), RIGHT, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            GrowArrow(v_vector),
            FadeIn(v_vector_label, shift=SMALL_BUFF*RIGHT),
            lag_ratio=0.4, run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        dist_brace_label = MathTex('v', 'T_B').scale(0.64).next_to(ground, DOWN, buff=SMALL_BUFF)
        dist_brace_label[0].set(color=VELOCITY_COLOR)
        dist_brace_label[1].set(color=TIME_COLOR)

        self.play(AnimationGroup(
            FadeOut(v_vector),
            v_vector_label.animate.move_to(dist_brace_label[0]),
            FadeOut(stopwatch, shift=SMALL_BUFF*DOWN),
            time_B_label.animate.move_to(dist_brace_label[1]),
            lag_ratio=0.05, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(v_vector_label, time_B_label)
        self.add(*dist_brace_label)

        self.wait(1.5)

        h_brace = BraceBetweenPoints(ORIGIN, 1.25*brace_height*DOWN).scale(0.8)
        h_brace.next_to(path_trace_line, LEFT, buff=SMALL_BUFF).shift(light.height/2*DOWN)
        h_brace_label = MathTex('h').scale(0.64).next_to(h_brace, LEFT, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            GrowFromEdge(h_brace, UP),
            FadeIn(h_brace_label, target_position=h_brace.get_left()),
            lag_ratio=0.5, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        path_brace = BraceBetweenPoints(ORIGIN, 1.25*path_trace_line.get_length()*LEFT).scale(0.8)
        path_angle = -np.arctan(h/(v*T_B))
        path_brace.rotate(path_angle).move_to(path_trace_line)
        path_brace.shift(Circle(0.15).point_at_angle(path_angle + PI/2))

        path_brace_label = MathTex(r'\sqrt{', 'h', '^2', '+', '(', 'v', 'T_B', ')', '^2}').scale(0.64)
        path_brace_label[5].set(color=VELOCITY_COLOR)
        path_brace_label[6].set(color=TIME_COLOR)
        path_brace_label.move_to(path_brace.get_center(), aligned_edge=LEFT).shift(0.432*UP)

        self.play(AnimationGroup(
            GrowFromCenter(path_brace),
            AnimationGroup(
                FadeIn(path_brace_label[0], shift=0),
                h_brace_label.animate.move_to(path_brace_label[1]),
                AnimationGroup(
                    dist_brace_label[0].animate.move_to(path_brace_label[5]),
                    dist_brace_label[1].animate.move_to(path_brace_label[6]),
                ),
                *[FadeIn(path_brace_label[i], shift=SMALL_BUFF*UP) for i in [2, 3, 4, 7, 8]],
                FadeOut(dist_brace),
                FadeOut(h_brace),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.3, run_time=5.5
        ))

        self.remove(*dist_brace_label, *h_brace_label)
        self.add(*path_brace_label)

        self.wait()

        c_vector = Vector(
            1/3 * (path_trace_line.get_end() - path_trace_line.get_start()),
            color=LIGHT_COLOR, stroke_width=4, max_tip_length_to_length_ratio=0.15,
        )
        c_vector.shift(light.get_center() + 4*light.height*DOWN)
        c_vector_label = MathTex('c', color=LIGHT_COLOR).scale(0.64).next_to(c_vector.get_tip(), DOWN, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            GrowArrow(c_vector),
            FadeIn(c_vector_label, shift=SMALL_BUFF*DOWN),
            lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        time_B_lhs_eq = MathTex('T_B', '=').scale(0.64)
        time_B_lhs_eq[0].set(color=TIME_COLOR)
        eq_buff = time_B_lhs_eq[1].get_left()[0] - time_B_lhs_eq[0].get_right()[0]
        time_B_eq_over = Line(ORIGIN, path_brace_label.width*RIGHT, stroke_width=1.5).next_to(time_B_lhs_eq, buff=eq_buff)
        Group(time_B_lhs_eq, time_B_eq_over).set_x(0).set_y(time_A_eq.get_y())

        self.play(AnimationGroup(
            FadeIn(time_B_lhs_eq, shift=0.2*LEFT),
            path_brace_label.animate.next_to(time_B_eq_over, UP, buff=SMALL_BUFF),
            c_vector_label.animate.next_to(time_B_eq_over, DOWN, buff=SMALL_BUFF),
            FadeIn(time_B_eq_over, shift=0.2*RIGHT),
            lag_ratio=0.05, run_time=3.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            c_vector_label.animate.next_to(time_B_lhs_eq, LEFT, buff=0).shift(0.01*DOWN + 0.02*LEFT),
            FadeOut(time_B_eq_over),
            path_brace_label.animate.set_y(time_B_lhs_eq.get_y()),
            FadeOut(path_brace),
            FadeOut(c_vector),
            lag_ratio=0.05, run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        time_B_eq_1 = MathTex('(', 'c', 'T_B', ')', '^2', '=', 'h', '^2', '+', '(', 'v', 'T_B', ')', '^2').scale(0.64)
        time_B_eq_1.set_x(0).set_y(time_A_eq.get_y())
        time_B_eq_1[1].set(color=LIGHT_COLOR)
        time_B_eq_1[2].set(color=TIME_COLOR)
        time_B_eq_1[10].set(color=VELOCITY_COLOR)
        time_B_eq_1[11].set(color=TIME_COLOR)

        self.play(AnimationGroup(
            c_vector_label.animate.move_to(time_B_eq_1[1]),
            time_B_lhs_eq[0].animate.move_to(time_B_eq_1[2]),
            time_B_lhs_eq[1].animate.move_to(time_B_eq_1[5]),
            *[path_brace_label[i].animate.move_to(time_B_eq_1[i+5]) for i in range(1, 9)],
            *[FadeIn(time_B_eq_1[i]) for i in [0, 3, 4, 5]],
            FadeOut(path_brace_label[0]),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(c_vector_label, *path_brace_label, *time_B_lhs_eq)
        self.add(*time_B_eq_1)

        self.wait(0.8)

        time_B_eq_2 = MathTex('c', '^2', 'T', '^2', '_B', '=', 'h', '^2', '+', 'v', '^2', 'T', '^2', '_B').scale(0.64)
        time_B_eq_2.move_to(time_B_eq_1)
        time_B_eq_2[0].set(color=LIGHT_COLOR)
        time_B_eq_2[2:5:2].set(color=TIME_COLOR)
        time_B_eq_2[9].set(color=VELOCITY_COLOR)
        time_B_eq_2[11:14:2].set(color=TIME_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_1, time_B_eq_2),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        time_B_eq_3 = MathTex('c', '^2', 'T', '^2', '_B', '-', 'v', '^2', 'T', '^2', '_B', '=', 'h', '^2').scale(0.64)
        time_B_eq_3.move_to(time_B_eq_2)
        time_B_eq_3[0].set(color=LIGHT_COLOR)
        time_B_eq_3[2:5:2].set(color=TIME_COLOR)
        time_B_eq_3[6].set(color=VELOCITY_COLOR)
        time_B_eq_3[8:11:2].set(color=TIME_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_2, time_B_eq_3),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        time_B_eq_4 = MathTex('(', 'c', '^2', '-', 'v', '^2', ')', 'T', '^2', '_B', '=', 'h', '^2').scale(0.64)
        time_B_eq_4.move_to(time_B_eq_3)
        time_B_eq_4[1].set(color=LIGHT_COLOR)
        time_B_eq_4[4].set(color=VELOCITY_COLOR)
        time_B_eq_4[7:10:2].set(color=TIME_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_3, time_B_eq_4),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        time_B_eq_5 = MathTex(r'\left(', '1', '-', '{v', '^2', r'\over', 'c', '^2}', r'\right)', 'T', '^2', '_B', '=', '{h', '^2', r'\over', 'c', '^2}').scale(0.64)
        time_B_eq_5.move_to(time_B_eq_4)
        time_B_eq_5[3].set(color=VELOCITY_COLOR)
        time_B_eq_5[6].set(color=LIGHT_COLOR)
        time_B_eq_5[9:12:2].set(color=TIME_COLOR)
        time_B_eq_5[16].set(color=LIGHT_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_4, time_B_eq_5),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        time_B_eq_6 = MathTex('T', '^2', '_B', '=', '{1', r'\over', '1', '-', '{v', '^2', r'\over', 'c', '^2}}', '{h', '^2', r'\over', 'c', '^2}').scale(0.64)
        time_B_eq_6.move_to(time_B_eq_5)
        time_B_eq_6[0:3:2].set(color=TIME_COLOR)
        time_B_eq_6[8].set(color=VELOCITY_COLOR)
        time_B_eq_6[11].set(color=LIGHT_COLOR)
        time_B_eq_6[16].set(color=LIGHT_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_5, time_B_eq_6),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        time_B_eq_7 = MathTex('T_B', '=', '{1', r'\over', r'\sqrt{', '1', '-', '{v', '^2', r'\over', 'c', '^2}}}', '{h', r'\over', 'c}').scale(0.64)
        time_B_eq_7.move_to(time_B_eq_6)
        time_B_eq_7[0].set(color=TIME_COLOR)
        time_B_eq_7[7].set(color=VELOCITY_COLOR)
        time_B_eq_7[10].set(color=LIGHT_COLOR)
        time_B_eq_7[14].set(color=LIGHT_COLOR)

        self.play(
            TransformMatchingTex(time_B_eq_6, time_B_eq_7),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        self.play(AnimationGroup(
            time_A_eq.animate.scale(time_B_eq_7[0].height/time_A_eq[0].height).set_x(-0.5, RIGHT).set_y(time_B_eq_7[1].get_y()),
            time_B_eq_7.animate.set_x(0.5, LEFT),
            run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            Circumscribe(time_A_eq[2:]),
            Circumscribe(time_B_eq_7[12:]),
            run_time=1.5
        ))

        self.wait(1.5)

        time_B_eq_8 = MathTex('T_B', '=', '{1', r'\over', r'\sqrt{', '1', '-', '{v', '^2', r'\over', 'c', '^2}}}', 'T_A').scale(0.64)
        time_B_eq_8.move_to(time_B_eq_6)
        time_B_eq_8[0].set(color=TIME_COLOR)
        time_B_eq_8[7].set(color=VELOCITY_COLOR)
        time_B_eq_8[10].set(color=LIGHT_COLOR)
        time_B_eq_8[12].set(color=TIME_COLOR)

        self.play(AnimationGroup(
            TransformMatchingShapes(time_B_eq_7[:12], time_B_eq_8[:12]),
            TransformMatchingShapes(time_A_eq[0], time_B_eq_8[12]),
            FadeOut(time_B_eq_7[12:]),
            FadeOut(time_A_eq[1:]),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        greater_than = MathTex('>', 'T_A').scale(0.64).next_to(time_B_eq_8[12], buff=0.2)
        greater_than[1].set(color=TIME_COLOR)

        self.play(
            FadeIn(greater_than, shift=0.2*RIGHT),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(3)
