from manim import *

###########################

POSITION_COLOR = BLUE
VELOCITY_COLOR = RED
ACCELERATION_COLOR = YELLOW

###########################

class TimePlot(Group):
    def __init__(self, func,
                 x_range, y_range,
                 x_length, y_length,
                 func_name=None,
                 axes_config={}, dot_config={}, curve_config={}, **kwargs):
        super().__init__(**kwargs)
        self._axes_config = axes_config
        self._dot_config = dot_config
        self._curve_config = curve_config
        self._kwargs = kwargs

        self.axes = Axes(
            x_range=x_range,
            y_range=y_range,
            x_length=x_length,
            y_length=y_length,
            axis_config={
                'include_ticks': False,
                'tip_width': 0.1,
                'tip_height': 0.1,
            },
            **axes_config,
            **kwargs
        )

        self.x_axis_label = self.axes.get_x_axis_label('t', edge=RIGHT, direction=RIGHT).scale(0.5)
        self.axes.add(self.x_axis_label)

        if func_name is not None and func_name != '':
            self.y_axis_label = self.axes.get_y_axis_label(func_name, edge=UP, direction=UP)
            self.y_axis_label.set(**kwargs).scale(0.5)
            self.axes.add(self.y_axis_label)

        self.dot = Dot(radius=0.05, z_index=100, **dot_config, **kwargs)
        self.dot.move_to(self.axes.c2p(0, func(0), 0))

        self.func = func
        self._hide_curve = False
        self._clip_curve = False
        self.curve = self.axes.plot(self.func, z_index=50, **curve_config, **kwargs)

        self.add(self.axes, self.dot, self.curve)

    def move_dot(self, t):
        self.dot.move_to(self.axes.c2p(t, self.func(t), 0))
        if self._clip_curve == True:
            self.curve.become(
                self.axes.plot(self.func, x_range=[0, t], **self._curve_config, **self._kwargs)
            )

    def hide_curve(self):
        self._hide_curve = True
        self._clip_curve = False
        self.curve.set(stroke_width=0)
        

    def clip_curve(self):
        self._hide_curve = False
        self._clip_curve = True
        t = self.axes.p2c(self.dot.get_center())[0]
        self.curve.become(
            self.axes.plot(self.func, x_range=[0, t], **self._curve_config, **self._kwargs)
        )

    def show_curve(self):
        self._hide_curve = False
        self._clip_curve = False
        self.curve.become(
            self.axes.plot(self.func, **self._curve_config, **self._kwargs)
        )
        

class ProjectileMotionScene(Scene):
    def construct(self):
        time = ValueTracker(0)

        x_0, y_0 = 0, 2
        g = 9.82

        d, h = 7.24, 3.05

        v_0_x = 5
        tau = d / v_0_x
        v_0_y = (h - y_0 + 1/2 * g * tau**2) / tau
        
        y_max = (v_0_y**2 + 2*g*y_0) / (2*g)

        def x(t):
            return x_0 + v_0_x * t

        def y(t):
            return y_0 + v_0_y * t - 1/2 * g * t**2

        def v_x(t):
            return v_0_x

        def v_y(t):
            return v_0_y - g * t

        def a_x(t):
            # Curve on plot does not render properly if this is 0 for some reason
            return 1e-6

        def a_y(t):
            return -g

        fill_conf = {'fill_color': WHITE, 'fill_opacity': 1, 'stroke_width': 0}
        ind_width = 0.2
        indicators = {
            'play': Triangle(**fill_conf).rotate(30*DEGREES).scale_to_fit_width(ind_width),
            'pause': VGroup(
                Rectangle(height=2*np.sqrt(3), width=1, **fill_conf),
                Rectangle(height=2*np.sqrt(3), width=1, **fill_conf),
            ).arrange(buff=1).scale_to_fit_width(ind_width),
            'rewind': VGroup(
                Triangle(**fill_conf).rotate(-30*DEGREES).scale_to_fit_width(ind_width),
                Triangle(**fill_conf).rotate(-30*DEGREES).scale_to_fit_width(ind_width),
            ).arrange(buff=0),
            'speed': VGroup(
                Triangle(**fill_conf).rotate(30*DEGREES).scale_to_fit_width(ind_width),
                Triangle(**fill_conf).rotate(30*DEGREES).scale_to_fit_width(ind_width),
            ).arrange(buff=0)
        }

        for ind in indicators.values():
            ind.to_corner(UL, buff=MED_LARGE_BUFF)

        cur_indicator = indicators['pause'].copy()

        court = Axes(
            x_range=[x_0, d],
            y_range=[0, y_max],
            x_length=d - x_0,
            y_length=y_max,
            y_axis_config={
                'stroke_width': 0
            },
            axis_config={
                'include_ticks': False,
            },
            tips=False
        )

        hoop = Line(ORIGIN, 0.4572*RIGHT, z_index=200).move_to(court.c2p(d, h, 0))
        brace_top = Line(ORIGIN, 0.1524*RIGHT).next_to(hoop, buff=0)
        brace = Polygon(
            ORIGIN, RIGHT, RIGHT+DOWN, 1/3*RIGHT+DOWN,
            stroke_width=0, fill_color=WHITE, fill_opacity=1
        ).scale_to_fit_width(0.1524).move_to(brace_top, aligned_edge=UP)
        backboard = Rectangle(
            height=1.0668, width=0.0508,
            stroke_width=0, fill_color=WHITE, fill_opacity=1
        ).move_to(court.c2p(d + 0.381, h - 0.3048, 0), aligned_edge=DL)
        net = Polygon(
            ORIGIN, RIGHT, 0.8*RIGHT+DOWN, 0.2*RIGHT+DOWN,
            stroke_width=2, color=WHITE
        ).scale_to_fit_width(0.4572).move_to(hoop, aligned_edge=UP)
        basket = Group(hoop, backboard, brace_top, net)

        ball = Dot(radius=0.12, z_index=100).move_to(court.c2p(x_0, y_0, 0))

        court_and_co = Group(court, basket)

        def ball_updater(mob):
            t = time.get_value()
            mob.move_to(court.c2p(x(t), y(t), 0))

        ball.add_updater(ball_updater)

        self.add(court_and_co, ball)
        self.add(cur_indicator)

        vector_config = {
            'tip_length': 0.2,
            'max_tip_length_to_length_ratio': 1,
            'stroke_width': 5,
            'max_stroke_width_to_length_ratio': 1000
        }

        v_0_vector = Vector(
            0.15*(v_0_x*RIGHT + v_0_y*UP),
            color=VELOCITY_COLOR, **vector_config
        ).shift(ball.get_center())
        v_0_vector_label = MathTex('v_0', font_size=30, color=VELOCITY_COLOR)
        v_0_vector_label.next_to(v_0_vector.get_tip(), UP)

        self.play(AnimationGroup(
            GrowArrow(v_0_vector),
            FadeIn(v_0_vector_label, target_position=v_0_vector.get_tip()),
            lag_ratio=0.5, run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))
        
        self.wait()
        
        cur_indicator.become(indicators['play'])
        self.play(AnimationGroup(
            FadeOut(v_0_vector, v_0_vector_label, run_time=0.5*tau, rate_func=rate_functions.ease_out_cubic),
            time.animate(run_time=tau, rate_func=rate_functions.linear).set_value(tau),
        ))
        cur_indicator.become(indicators['pause'])

        self.wait(2)

        x_line = DashedLine(
            court.c2p(0, h, 0), court.c2p(d, h, 0),
            dash_length=0.06, dashed_ratio=0.4, stroke_width=3, color=POSITION_COLOR
        )
        y_line = DashedLine(
            court.c2p(d, 0, 0), court.c2p(d, y_max, 0),
            dash_length=0.06, dashed_ratio=0.4, stroke_width=3, color=POSITION_COLOR
        )

        x_line_label = MathTex('x', font_size=30, color=POSITION_COLOR)
        y_line_label = MathTex('y', font_size=30, color=POSITION_COLOR)

        x_line_label_width = x_line_label.width

        def x_line_updater(mob):
            t = time.get_value()
            try:
                mob.put_start_and_end_on(court.c2p(0, y(t), 0), court.c2p(d, y(t), 0))
            except:
                pass
            for line in mob.submobjects:
                if line.get_x() > ball.get_x():
                    line.set(color=BLACK)
                else:
                    line.set(color=POSITION_COLOR)

        def y_line_updater(mob):
            t = time.get_value()
            try:
                mob.put_start_and_end_on(court.c2p(x(t), 0, 0), court.c2p(x(t), y_max, 0))
            except:
                pass
            for line in mob.submobjects:
                if line.get_y() > ball.get_y():
                    line.set(color=BLACK)
                else:
                    line.set(color=POSITION_COLOR)

        def x_line_label_updater(mob):
            t = time.get_value()
            if x(t) == 0:
                mob.move_to([999, 999, 0])
            else:
                mob.scale_to_fit_width(
                    min(x_line_label_width, 0.3*x(t))
                )
                mob.move_to(court.c2p(0.5*x(t), y(t), 0) + 0.15*UP, aligned_edge=DOWN)

        def y_line_label_updater(mob):
            t = time.get_value()
            mob.move_to(
                court.c2p(x(t), 0.5*y(t), 0)
                    + 0.15*RIGHT, aligned_edge=LEFT
            )

        x_line_updater(x_line)
        x_line_label_updater(x_line_label)
        y_line_updater(y_line)
        y_line_label_updater(y_line_label)

        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromPoint(x_line, ball),
                FadeIn(x_line_label, shift=0.15*UP),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(y_line, ball),
                FadeIn(y_line_label, shift=0.15*RIGHT),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.1, run_time=4
        ))

        x_line.add_updater(x_line_updater)
        x_line_label.add_updater(x_line_label_updater)
        y_line.add_updater(y_line_updater)
        y_line_label.add_updater(y_line_label_updater)

        plot_config = {
            'x_range': [0, tau],
            'x_length': 1.2,
            'y_length': 1.6,
        }

        s_plots = Group(
            TimePlot(
                x, func_name='x',
                y_range=[0, 1.2*d],
                color=POSITION_COLOR, **plot_config
            ),
            TimePlot(
                y, func_name='y',
                y_range=[0, 1.2*d],
                color=POSITION_COLOR, **plot_config
            )
        )

        def prepare_plots(plots):
            for plot in plots:
                plot.axes.get_y_axis().add_labels({0: '$0$'}, font_size=15, buff=0.15)
                plot.add_updater(lambda mob: mob.move_dot(time.get_value()))
                plot.hide_curve()

        prepare_plots(s_plots)

        s_plots.scale(1.5).arrange(DOWN, buff=MED_LARGE_BUFF).next_to(court_and_co, buff=LARGE_BUFF)
        
        on_screen = Group(court_and_co, s_plots)
        shift_val = on_screen.copy().center().get_center() - on_screen.get_center()

        s_plots.shift(shift_val)

        self.play(AnimationGroup(
            court_and_co.animate.shift(shift_val),
            AnimationGroup(FadeIn(s_plots[0], shift=shift_val)),
            AnimationGroup(FadeIn(s_plots[1], shift=shift_val)),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        cur_indicator.become(indicators['rewind'])
        self.play(time.animate.set_value(0.3*tau), run_time=3)
        cur_indicator.become(indicators['speed'])
        self.play(time.animate.set_value(0.6*tau), run_time=3)
        cur_indicator.become(indicators['rewind'])
        self.play(time.animate.set_value(0), run_time=3)
        cur_indicator.become(indicators['pause'])

        self.wait(2)

        y_0_brace = BraceBetweenPoints(court.c2p(0, 2*y_0, 0), court.get_origin()).scale(0.5, about_edge=DR)
        y_0_brace_label = MathTex('y_0', font_size=30).next_to(y_0_brace, LEFT, buff=0.15)

        s_plots.scale(1/1.5)
        y_0_line = s_plots[1].axes.get_horizontal_line(s_plots[1].axes.c2p(tau, y_0, 0), color=GRAY)
        y_0_line_label = MathTex('y_0', font_size=15).next_to(y_0_line, LEFT, buff=0.15)
        Group(s_plots, y_0_line, y_0_line_label).scale(1.5, about_point=s_plots.get_center())

        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromCenter(y_0_brace),
                FadeIn(y_0_brace_label, target_position=y_0_brace),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(y_0_line, s_plots[1].dot),
                FadeIn(y_0_line_label, target_position=y_0_line.get_left()),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=3
        ))

        s_plots[1].add(y_0_line, y_0_line_label)

        self.wait()

        d_brace = BraceBetweenPoints(court.get_origin(), court.c2p(2*d, 0, 0)).scale(0.5, about_edge=UL)
        d_brace_label = MathTex('d', font_size=30).next_to(d_brace, DOWN, buff=0.15)

        s_plots.scale(1/1.5)
        d_line = s_plots[0].axes.get_horizontal_line(s_plots[0].axes.c2p(tau, d, 0), color=GRAY)
        d_line_label = MathTex('d', font_size=15).next_to(d_line, LEFT, buff=0.15)
        Group(s_plots, d_line, d_line_label).scale(1.5, about_point=s_plots.get_center())

        cur_indicator.become(indicators['speed'])
        self.play(time.animate.set_value(tau), run_time=2)
        cur_indicator.become(indicators['pause'])
        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromCenter(d_brace),
                FadeIn(d_brace_label, target_position=d_brace),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(d_line, s_plots[0].axes.c2p(tau, d, 0)),
                FadeIn(d_line_label, target_position=d_line.get_left()),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=3
        ))

        s_plots[0].add(d_line, d_line_label)

        self.wait()

        h_brace_x = court.p2c(brace_top.get_right())[0]
        h_brace = BraceBetweenPoints(court.c2p(h_brace_x, 0, 0), court.c2p(h_brace_x, 2*h, 0)).scale(0.5, about_edge=DL)
        h_brace_label = MathTex('h', font_size=30).next_to(h_brace, RIGHT, buff=0.15)

        s_plots.scale(1/1.5)
        h_line = s_plots[1].axes.get_horizontal_line(s_plots[1].axes.c2p(tau, h, 0), color=GRAY)
        h_line_label = MathTex('h', font_size=15).next_to(h_line, LEFT, buff=0.15)
        Group(s_plots, h_line, h_line_label).scale(1.5, about_point=s_plots.get_center())

        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromCenter(h_brace),
                FadeIn(h_brace_label, target_position=h_brace),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(h_line, s_plots[1].dot),
                FadeIn(h_line_label, target_position=h_line.get_left()),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=3
        ))
        
        s_plots[1].add(h_line, h_line_label)

        self.wait()

        cur_indicator.become(indicators['rewind'])
        self.play(time.animate.set_value(0), run_time=2)
        cur_indicator.become(indicators['pause'])

        self.wait()

        for plot in s_plots:
            plot.clip_curve()

        cur_indicator.become(indicators['play'])
        ind_speed = Tex('0.5x', font_size=20).next_to(cur_indicator)
        self.play(
            AnimationGroup(
                FadeIn(ind_speed, target_position=cur_indicator.get_right()),
                run_time=0.75, rate_func=rate_functions.ease_out_cubic
            ),
            time.animate(run_time=2*tau, rate_func=rate_functions.linear).set_value(tau)
        )
        cur_indicator.become(indicators['pause'])
        self.play(
            FadeOut(ind_speed, target_position=cur_indicator.get_right()),
            run_time=0.75, rate_func=rate_functions.ease_out_cubic
        )

        for plot in s_plots:
            plot.show_curve()

        self.wait()

        shift_val = -shift_val

        self.play(
            court_and_co.animate.shift(shift_val),
            FadeOut(
                x_line, x_line_label,
                y_line, y_line_label,
                y_0_brace, y_0_brace_label,
                d_brace, d_brace_label,
                h_brace, h_brace_label,
                s_plots,
                shift=shift_val
            ),
            time.animate.set_value(0.2*tau),
            run_time=4
        )

        self.wait(2.5)

        v_vector = Vector()
        v_vector_label = MathTex('v', font_size=30, color=VELOCITY_COLOR, z_index=50)

        def v_vector_updater(mob):
            t = time.get_value()
            mob.become(
                Vector(
                    0.15*(v_x(t)*RIGHT + v_y(t)*UP),
                    color=VELOCITY_COLOR, z_index=50, **vector_config
                ).shift(ball.get_center())
            )

        def v_vector_label_updater(mob):
            mob.next_to(v_vector.get_tip(), UP)

        v_vector_updater(v_vector)
        v_vector_label_updater(v_vector_label)

        self.play(AnimationGroup(
            GrowArrow(v_vector),
            FadeIn(v_vector_label, target_position=v_vector.get_tip()),
            lag_ratio=0.5, run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        v_vector.add_updater(v_vector_updater)
        v_vector_label.add_updater(v_vector_label_updater)

        self.wait()

        v_vector_eq = MathTex('{{ v }} =', font_size=30)
        v_vector_eq.shift(v_vector_label.get_center() - v_vector_eq.submobjects[0].get_center())
        v_vector_eq_buff = v_vector_eq.submobjects[1].get_left()[0] - v_vector_eq.submobjects[0].get_right()[0]
        v_vector_matrix = Matrix(
            [['v_x'], ['v_y']],
            left_bracket='(', right_bracket=')',
            v_buff=0.5, bracket_v_buff=SMALL_BUFF, bracket_h_buff=SMALL_BUFF,
            stretch_brackets=False,
            element_to_mobject_config={'font_size': 30}
        ).set_column_colors(VELOCITY_COLOR).next_to(v_vector_eq, buff=v_vector_eq_buff)
        for bracket in v_vector_matrix.get_brackets():
            bracket.scale_to_fit_height(v_vector_matrix.get_entries().height + 2*SMALL_BUFF)

        self.play(
            FadeIn(
                v_vector_eq.submobjects[1], v_vector_matrix,
                shift=v_vector_eq_buff*RIGHT
            ),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        v_x_line = DashedLine(
            v_vector.get_start(), v_vector.get_corner(DR),
            stroke_width=2, color=VELOCITY_COLOR, z_index=50
        )
        v_y_line = Group(
            *DashedLine(
                ORIGIN, 0.15*v_0_y*UP,
                stroke_width=2, color=VELOCITY_COLOR
            ),
            *DashedLine(
                ORIGIN, 0.15*v_0_y*DOWN,
                stroke_width=2, color=VELOCITY_COLOR
            ),
            z_index=50
        )

        v_x_line_label = MathTex('v_x', font_size=20, color=VELOCITY_COLOR, z_index=50)
        v_y_line_label = MathTex('v_y', font_size=20, color=VELOCITY_COLOR, z_index=50)

        def v_x_line_updater(mob):
            mob.move_to(v_vector.get_start(), aligned_edge=LEFT)

        def v_y_line_updater(mob):
            mob.move_to([v_vector.get_end()[0], v_vector.get_start()[1], 0])
            for line in mob.submobjects:
                line_y = line.get_y() - v_vector.get_start()[1]
                end_y = v_vector.get_end()[1] - v_vector.get_start()[1]
                if line_y * end_y < 0 or abs(line_y) > abs(end_y):
                    line.set(color=BLACK)
                else:
                    line.set(color=VELOCITY_COLOR)

        def v_x_line_label_updater(mob):
            t = time.get_value()
            buff = 2 * rate_functions.ease_in_out_cubic((sorted([-1, v_y(t), 1])[1] + 1) / 2) - 1
            buff *= 0.15 + v_x_line_label.height/2
            if v_y(t) > 0:
                mob.move_to(v_vector.get_bottom()).shift(buff*DOWN)
            else:
                mob.move_to(v_vector.get_top()).shift(buff*DOWN)
            shift_x = mob.get_right()[0] - (backboard.get_left()[0] - SMALL_BUFF)
            if shift_x > 0:
                mob.shift(shift_x*LEFT)

        def v_y_line_label_updater(mob):
            mob.next_to(v_vector, RIGHT, buff=0.15)

        v_x_line_updater(v_x_line)
        v_y_line_updater(v_y_line)
        v_x_line_label_updater(v_x_line_label)
        v_y_line_label_updater(v_y_line_label)

        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromPoint(v_x_line, v_vector.get_start()),
                Transform(
                    v_vector_matrix.get_entries()[0], v_x_line_label,
                    replace_mobject_with_target_in_scene=True
                ),
                lag_ratio=0.1, run_time=4, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(v_y_line, [v_vector.get_end()[0], v_vector.get_start()[1], 0]),
                Transform(
                    v_vector_matrix.get_entries()[1], v_y_line_label,
                    replace_mobject_with_target_in_scene=True
                ),
                lag_ratio=0.1, run_time=4, rate_func=rate_functions.ease_out_cubic
            ),
            FadeOut(
                *v_vector_eq.submobjects[1:], v_vector_matrix.get_brackets(),
                run_time=1.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2
        ))

        v_x_line.add_updater(v_x_line_updater)
        v_y_line.add_updater(v_y_line_updater)
        v_x_line_label.add_updater(v_x_line_label_updater)
        v_y_line_label.add_updater(v_y_line_label_updater)

        self.wait()

        v_bound = max(abs(v_0_y), abs(v_0_y - g*tau))
        v_plots = Group(
            TimePlot(
                v_x, func_name='v_x',
                y_range=[-1.3*v_bound, 1.3*v_bound],
                color=VELOCITY_COLOR, **plot_config
            ),
            TimePlot(
                v_y, func_name='v_y',
                y_range=[-1.3*v_bound, 1.3*v_bound],
                color=VELOCITY_COLOR, **plot_config
            )
        )

        prepare_plots(v_plots)

        v_plots.scale(1.5).arrange(DOWN, buff=MED_LARGE_BUFF).next_to(court_and_co, buff=LARGE_BUFF)
        
        on_screen = Group(court_and_co, v_plots)
        shift_val = on_screen.copy().center().get_center() - on_screen.get_center()

        v_plots.shift(shift_val)

        self.play(AnimationGroup(
            court_and_co.animate.shift(shift_val),
            AnimationGroup(FadeIn(v_plots[0], shift=shift_val)),
            AnimationGroup(FadeIn(v_plots[1], shift=shift_val)),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        t_is_zero = MathTex('t = 0', font_size=30).next_to(indicators['rewind'])

        cur_indicator.become(indicators['rewind'])
        self.play(AnimationGroup(
            time.animate.set_value(0),
            FadeIn(t_is_zero, shift=MED_LARGE_BUFF*RIGHT),
            lag_ratio=0.3, run_time=1.5
        ))
        cur_indicator.become(indicators['pause'])
        
        self.wait(0.5)

        v_vector_eq = MathTex('{{ v }} = {{ v_0 }}', font_size=30)
        v_vector_eq.submobjects[2].set(color=VELOCITY_COLOR)
        v_vector_eq.shift(v_vector_label.get_center() - v_vector_eq.submobjects[0].get_center())

        v_0_x_eq = MathTex('{{ v_x }} = {{ v_{0,x} }}', font_size=20)
        v_0_x_eq.submobjects[2].set(color=VELOCITY_COLOR)
        v_0_x_eq_buff = v_0_x_eq.submobjects[1].get_left()[0] - v_0_x_eq.submobjects[0].get_right()[0]
        v_0_x_eq.shift(v_x_line_label.get_center() - v_0_x_eq.submobjects[0].get_center())

        v_0_y_eq = MathTex('{{ v_y }} = {{ v_{0,y} }}', font_size=20)
        v_0_y_eq.submobjects[2].set(color=VELOCITY_COLOR)
        v_0_y_eq_buff = v_0_y_eq.submobjects[1].get_left()[0] - v_0_y_eq.submobjects[0].get_right()[0]
        v_0_y_eq.shift(v_y_line_label.get_center() - v_0_y_eq.submobjects[0].get_center())

        v_plots.scale(1/1.5)
        v_0_x_line = v_plots[0].axes.get_horizontal_line(v_plots[0].axes.c2p(tau, v_0_x, 0), color=GRAY)
        v_0_x_line_label = MathTex('v_{0,x}', font_size=15).next_to(v_0_x_line, LEFT, buff=0.15)
        Group(v_plots, v_0_x_line, v_0_x_line_label).scale(1.5, about_point=v_plots.get_center())

        v_plots.scale(1/1.5)
        v_0_y_line = v_plots[1].axes.get_horizontal_line(v_plots[1].axes.c2p(tau, v_0_y, 0), color=GRAY)
        v_0_y_line_label = MathTex('v_{0,y}', font_size=15).next_to(v_0_y_line, LEFT, buff=0.15)
        Group(v_plots, v_0_y_line, v_0_y_line_label).scale(1.5, about_point=v_plots.get_center())

        self.play(AnimationGroup(
            AnimationGroup(
                FadeIn(*v_vector_eq.submobjects[1:], shift=v_vector_eq_buff*RIGHT),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(*v_0_x_eq.submobjects[1:], shift=v_0_x_eq_buff*RIGHT),
                AnimationGroup(
                    GrowFromPoint(v_0_x_line, v_plots[0].dot),
                    FadeIn(v_0_x_line_label, target_position=v_0_x_line.get_left()),
                    lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
                ),
                FadeIn(*v_0_y_eq.submobjects[1:], shift=v_0_y_eq_buff*RIGHT),
                AnimationGroup(
                    GrowFromPoint(v_0_y_line, v_plots[1].dot),
                    FadeIn(v_0_y_line_label, target_position=v_0_y_line.get_left()),
                    lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
                ),
                lag_ratio=0.3, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.6, run_time=8
        ))

        v_plots[0].add(v_0_x_line, v_0_x_line_label)
        v_plots[1].add(v_0_y_line, v_0_y_line_label)

        self.wait()
        self.play(
            FadeOut(
                t_is_zero, v_vector_label, 
                *v_vector_eq.submobjects[1:],
                *v_0_x_eq.submobjects[1:],
                *v_0_y_eq.submobjects[1:]
            )
        )
        
        for plot in v_plots:
            plot.clip_curve()

        cur_indicator.become(indicators['play'])
        self.play(
            AnimationGroup(
                FadeIn(ind_speed, target_position=cur_indicator.get_right()),
                run_time=0.75, rate_func=rate_functions.ease_out_cubic
            ),
            time.animate(run_time=2*tau, rate_func=rate_functions.linear).set_value(tau)
        )
        cur_indicator.become(indicators['pause'])
        self.play(
            FadeOut(ind_speed, target_position=cur_indicator.get_right()),
            run_time=0.75, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        for plot in v_plots:
            plot.show_curve()

        cur_indicator.become(indicators['rewind'])
        self.play(time.animate.set_value(0), run_time=2)
        cur_indicator.become(indicators['pause'])

        a_vector = Vector()
        a_vector_label = MathTex('a', font_size=30, color=ACCELERATION_COLOR)

        def a_vector_updater(mob):
            t = time.get_value()
            mob.become(
                Vector(
                    0.1*(a_x(t)*RIGHT + a_y(t)*UP),
                    color=ACCELERATION_COLOR, **vector_config
                ).shift(ball.get_center())
            )

        def a_vector_label_updater(mob):
            mob.next_to(a_vector.get_tip(), RIGHT)

        a_vector_updater(a_vector)
        a_vector_label_updater(a_vector_label)

        a_plots = Group(
            TimePlot(
                a_x, func_name='a_x',
                y_range=[-15, 15],
                color=ACCELERATION_COLOR, **plot_config
            ),
            TimePlot(
                a_y, func_name='a_y',
                y_range=[-15, 15],
                color=ACCELERATION_COLOR, **plot_config
            )
        )

        prepare_plots(a_plots)

        for plot in a_plots:
            plot.clip_curve()

        a_plots.scale(1.5).arrange(DOWN, buff=MED_LARGE_BUFF).next_to(court_and_co, buff=LARGE_BUFF)

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(
                    v_vector, v_vector_label,
                    v_x_line, v_x_line_label,
                    v_y_line, v_y_line_label
                ),
                run_time=1.5, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeOut(v_plots, shift=2*LEFT),
                run_time=3
            ),
            AnimationGroup(
                FadeIn(a_plots, shift=2*LEFT),
                run_time=3
            ),
            AnimationGroup(
                GrowArrow(a_vector),
                FadeIn(a_vector_label, target_position=a_vector.get_tip()),
                lag_ratio=0.5, run_time=4.5, rate_func=rate_functions.ease_out_cubic
            )
        ))

        self.wait(1.5)

        a_vector.add_updater(a_vector_updater)
        a_vector_label.add_updater(a_vector_label_updater)

        a_vector_eq = MathTex('{{ a }} =', font_size=30)
        a_vector_eq.shift(a_vector_label.get_center() - a_vector_eq.submobjects[0].get_center())
        a_vector_eq_buff = a_vector_eq.submobjects[1].get_left()[0] - a_vector_eq.submobjects[0].get_right()[0]
        a_vector_matrix = Matrix(
            [['0'], ['-g']],
            left_bracket='(', right_bracket=')',
            v_buff=0.5, bracket_v_buff=SMALL_BUFF, bracket_h_buff=SMALL_BUFF,
            stretch_brackets=False,
            element_to_mobject_config={'font_size': 30}
        ).set_column_colors(ACCELERATION_COLOR).next_to(a_vector_eq, buff=a_vector_eq_buff)

        for bracket in v_vector_matrix.get_brackets():
            bracket.scale_to_fit_height(v_vector_matrix.get_entries().height + 2*SMALL_BUFF)

        a_plots.scale(1/1.5)
        g_line = a_plots[1].axes.get_horizontal_line(a_plots[1].axes.c2p(tau, -g, 0), color=GRAY)
        g_line_label = MathTex('-g', font_size=15).next_to(g_line, LEFT, buff=0.15)
        Group(a_plots, g_line, g_line_label).scale(1.5, about_point=a_plots.get_center())
            
        self.play(AnimationGroup(
            FadeIn(
                a_vector_eq.submobjects[1], a_vector_matrix,
                shift=a_vector_eq_buff*RIGHT, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                GrowFromPoint(g_line, a_plots[1].dot),
                FadeIn(g_line_label, target_position=g_line.get_left()),
                lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.6, run_time=4.5
        ))

        a_plots[1].add(g_line, g_line_label)

        self.wait(0.75)

        self.play(
            FadeOut(a_vector_eq.submobjects[1], a_vector_matrix, shift=a_vector_eq_buff*LEFT),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        cur_indicator.become(indicators['play'])
        self.play(
            AnimationGroup(
                FadeIn(ind_speed, target_position=cur_indicator.get_right()),
                run_time=0.75, rate_func=rate_functions.ease_out_cubic
            ),
            time.animate(run_time=2*tau, rate_func=rate_functions.linear).set_value(tau)
        )
        cur_indicator.become(indicators['pause'])
        self.play(
            FadeOut(ind_speed, target_position=cur_indicator.get_right()),
            run_time=0.75, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(3)
        
        for plot in a_plots:
            plot.show_curve()

        cur_indicator.become(indicators['rewind'])
        self.play(time.animate.set_value(0), run_time=2)
        cur_indicator.become(indicators['pause'])

        self.wait(0.3)

        target_plots = Group(*s_plots, *v_plots, *a_plots).copy()
        target_plots.arrange_in_grid(rows=3, cols=2, buff=MED_LARGE_BUFF)
        target_plots.scale(0.75).next_to(court_and_co, buff=MED_LARGE_BUFF)

        on_screen = Group(court_and_co, target_plots)
        shift_val = on_screen.copy().center().get_center() - on_screen.get_center()
        target_plots.shift(shift_val)
        
        s_plots = target_plots[:2]
        v_plots = target_plots[2:4]
        plot_shift_val = (4 + MED_LARGE_BUFF - v_plots.get_bottom()[1])*UP
        Group(s_plots, v_plots).shift(plot_shift_val)

        self.play(AnimationGroup(
            AnimationGroup(
                AnimationGroup(
                    Transform(a_plots, target_plots[4:]),
                    run_time=3, rate_func=rate_functions.ease_out_cubic
                ),
                court_and_co.animate(run_time=3, rate_func=rate_functions.ease_out_cubic).shift(shift_val),
                v_plots.animate(run_time=3, rate_func=rate_functions.ease_out_cubic).shift(-plot_shift_val),
                s_plots.animate(run_time=3, rate_func=rate_functions.ease_out_cubic).shift(-plot_shift_val),
                lag_ratio=0.05
            ),
        ))

        mobs = [
            x_line, x_line_label,
            y_line, y_line_label,
            v_vector, 
            v_x_line, v_x_line_label,
            v_y_line, v_y_line_label,
        ]
        updaters = [
            x_line_updater, x_line_label_updater,
            y_line_updater, y_line_label_updater,
            v_vector_updater, 
            v_x_line_updater, v_x_line_label_updater,
            v_y_line_updater, v_y_line_label_updater,
        ]

        for mob, updater in zip(mobs, updaters):
            updater(mob)
            mob.remove_updater(updater)

        a_vector_label.remove_updater(a_vector_label_updater)

        self.play(AnimationGroup(
            AnimationGroup(
                AnimationGroup(
                    GrowArrow(v_vector),
                    run_time=1.5, rate_func=rate_functions.ease_out_cubic
                ),
                AnimationGroup(
                    AnimationGroup(
                        GrowFromPoint(v_x_line, v_vector.get_start()),
                        FadeIn(v_x_line_label, target_position=v_x_line),
                        lag_ratio=0.3, run_time=3, rate_func=rate_functions.ease_out_cubic
                    ),
                    AnimationGroup(
                        GrowFromPoint(v_y_line, [v_vector.get_end()[0], v_vector.get_start()[1], 0]),
                        FadeIn(v_y_line_label, target_position=v_vector.get_right()),
                        lag_ratio=0.3, run_time=3, rate_func=rate_functions.ease_out_cubic
                    ),
                    lag_ratio=0.2
                )
            ),

            a_vector_label.animate(rate_func=rate_functions.ease_out_cubic).next_to(a_vector.get_tip(), LEFT, buff=0.15),

            AnimationGroup(
                AnimationGroup(
                    GrowFromPoint(x_line, ball),
                    FadeIn(x_line_label, shift=0.15*UP),
                    lag_ratio=0.5, run_time=1.5, rate_func=rate_functions.ease_out_cubic
                ),
                AnimationGroup(
                    GrowFromPoint(y_line, ball),
                    FadeIn(y_line_label, shift=0.15*RIGHT),
                    lag_ratio=0.5, rate_func=rate_functions.ease_out_cubic
                ),
                lag_ratio=0.1, run_time=4
            ),
            lag_ratio=0.2
        ))

        def a_vector_label_updater(mob):
            mob.next_to(a_vector.get_tip(), LEFT, buff=0.15)

        for mob, updater in zip(mobs, updaters):
            mob.add_updater(updater)
        
        a_vector_label.add_updater(a_vector_label_updater)

        self.wait(1.5)

        cur_indicator.become(indicators['play'])
        self.play(
            AnimationGroup(
                FadeIn(ind_speed, target_position=cur_indicator.get_right()),
                run_time=0.75, rate_func=rate_functions.ease_out_cubic
            ),
            time.animate(run_time=2*tau, rate_func=rate_functions.linear).set_value(tau)
        )
        cur_indicator.become(indicators['pause'])
        self.play(
            FadeOut(ind_speed, target_position=cur_indicator.get_right()),
            run_time=0.75, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(3)
