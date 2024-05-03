from manim import *

#####################

POSITIVE_COLOR = RED
NEGATIVE_COLOR = BLUE
CHARGE_COLOR = YELLOW

#####################

class FieldLine(VMobject):
    def __init__(self, field_func, start, stop_cond=None, max_length=10, step_size=1e-2, **kwargs):
        super().__init__(**kwargs)
        self.start_new_path(start)
        pos = start
        for _ in np.arange(max_length/step_size):
            field_val = field_func(pos)
            if field_val is None:
                break
            pos += step_size * field_val / np.linalg.norm(field_val)
            if stop_cond is not None and stop_cond(pos):
                break
            self.add_line_to(pos)


class ElectricFieldLinesScene(Scene):
    def construct(self):
        negative = Group(Dot(radius=0.2, color=NEGATIVE_COLOR), MathTex('-', color=BLACK).scale(0.8)).shift(2*LEFT)
        positive = Group(Dot(radius=0.2, color=POSITIVE_COLOR), MathTex('+', color=BLACK).scale(0.8)).shift(2*RIGHT)

        self.add(positive, negative)

        charge = Group(Dot(radius=0.25, color=CHARGE_COLOR), MathTex('+', color=BLACK)).scale(0.5)
        charge.shift(UP)
        
        self.wait()

        self.play(
            FadeIn(charge, shift=0.15*UP),
            run_time=2.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        def neg_field(v):
            r_minus = negative.get_center() - v
            if np.linalg.norm(r_minus) == 0:
                return ORIGIN
            return (1 / np.linalg.norm(r_minus)**2) * (r_minus / np.linalg.norm(r_minus))

        vec_scale = 3
        vec_config = {
            'max_tip_length_to_length_ratio': 0.25,
            'max_stroke_width_to_length_ratio': 99,
            'stroke_width': 6,
        }

        def on_charge(angle):
            return charge.get_center() + Circle(0.1375).point_at_angle(angle)

        neg_vec = Vector(vec_scale*neg_field(charge.get_center()), color=NEGATIVE_COLOR, **vec_config).shift(on_charge(angle_of_vector(neg_field(charge.get_center()))))

        self.play(
            GrowArrow(neg_vec),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(0.5)

        force_eq = MathTex('F', '=', 'k_B', '{q', 'Q', r'\over', 'r', '^2}').scale(0.8).to_corner(UL, buff=LARGE_BUFF)
        force_eq[0:5:4].set(color=NEGATIVE_COLOR)
        force_eq[3].set(color=CHARGE_COLOR)

        self.play(
            FadeIn(force_eq, shift=LARGE_BUFF*RIGHT),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        def pos_field(v):
            r_plus = positive.get_center() - v
            if np.linalg.norm(r_plus) == 0:
                return ORIGIN
            return (1 / np.linalg.norm(r_plus)**2) * (-r_plus / np.linalg.norm(r_plus))

        pos_vec = Vector(vec_scale*pos_field(charge.get_center()), color=POSITIVE_COLOR, **vec_config).shift(on_charge(angle_of_vector(pos_field(charge.get_center()))))
        
        self.play(
            GrowArrow(pos_vec),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        self.play(
            force_eq[0:5:4].animate.set(color=POSITIVE_COLOR),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait(3)

        self.play(
            FadeOut(force_eq, shift=LARGE_BUFF*LEFT),
            run_time=1.5, rate_func=rate_functions.ease_out_cubic
        )

        def total_field(v):
            neg = neg_field(v)
            pos = pos_field(v)
            if np.linalg.norm(neg) == 0 or np.linalg.norm(pos) == 0:
                return ORIGIN
            return neg + pos

        force_vec = Vector(vec_scale*total_field(charge.get_center()), **vec_config).shift(on_charge(angle_of_vector(total_field(charge.get_center()))))

        self.play(
            GrowArrow(force_vec),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        neg_vec.add_updater(lambda mob: mob.put_start_and_end_on(ORIGIN, vec_scale*neg_field(charge.get_center())).shift(on_charge(angle_of_vector(neg_field(charge.get_center())))))
        pos_vec.add_updater(lambda mob: mob.put_start_and_end_on(ORIGIN, vec_scale*pos_field(charge.get_center())).shift(on_charge(angle_of_vector(pos_field(charge.get_center())))))
        force_vec.add_updater(lambda mob: mob.put_start_and_end_on(ORIGIN, vec_scale*total_field(charge.get_center())).shift(on_charge(angle_of_vector(total_field(charge.get_center())))))
        
        charge_path = VMobject()
        charge_path.append_vectorized_mobject(Line(charge.get_center(), ORIGIN))
        charge_path.append_vectorized_mobject(Arc(radius=2, start_angle=PI, angle=PI, arc_center=positive.get_center()))
        charge_path.append_vectorized_mobject(ParametricFunction(lambda t: 4*np.cos(t)*(RIGHT + np.sin(t)*UP), t_range=[0, 1.5*PI]))

        self.play(
            MoveAlongPath(charge, charge_path),
            run_time=16, rate_func=rate_functions.smooth
        )

        self.wait(2)

        self.play(
            FadeOut(Group(charge, neg_vec, pos_vec, force_vec), shift=0.15*UP),
            run_time=2, rate_func=rate_functions.ease_in_cubic
        )

        self.wait()

        vec_field_config = {
            'length_func': lambda _: 0.35,
            'min_color_scheme_value': 0,
            'max_color_scheme_value': 0.15,
        }

        neg_vec_field = ArrowVectorField(neg_field, colors=[DARK_GRAY, NEGATIVE_COLOR], **vec_field_config)
        self.play(FadeIn(neg_vec_field), run_time=5, rate_func=rate_functions.ease_out_cubic)

        self.wait(2)

        self.play(neg_vec_field.animate.set_opacity(0.33))

        pos_vec_field = ArrowVectorField(pos_field, colors=[DARK_GRAY, POSITIVE_COLOR], **vec_field_config)
        self.play(FadeIn(pos_vec_field), run_time=5, rate_func=rate_functions.ease_out_cubic)

        self.wait(2)

        self.play(neg_vec_field.animate.set_opacity(1))

        self.wait(0.5)

        total_vec_field = ArrowVectorField(total_field, colors=[DARK_GRAY, WHITE], **vec_field_config)
        self.play(AnimationGroup(
            *[ReplacementTransform(pos_arrow, total_arrow) for pos_arrow, total_arrow in zip(pos_vec_field, total_vec_field)],
            *[ReplacementTransform(neg_arrow, total_arrow) for neg_arrow, total_arrow in zip(neg_vec_field, total_vec_field)],
            run_time=6, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.wait(2)

        field_lines_full = []
        field_lines_partial = []

        for angle in np.linspace(0, 2*PI, 20, endpoint=False):
            source = positive.get_center() + Circle(0.22).point_at_angle(angle)
            field_line_from = FieldLine(total_field, source, lambda pos: pos[0] < 0, stroke_width=1)

            sink = negative.get_center() + Circle(0.22).point_at_angle(PI - angle)
            field_line_to = FieldLine(lambda v: -total_field(v), sink, lambda pos: pos[0] > 0, stroke_width=1)
            field_line_to.reverse_direction()

            if PI/2 <= angle <= 3*PI/2:
                field_lines_full.append((field_line_from, field_line_to))
            else:
                field_lines_partial.append((field_line_from, field_line_to))

        arrows_full = []

        for field_line_from, _ in field_lines_full:
            arrow = ArrowTriangleFilledTip(color=WHITE).scale(0.4).move_to(field_line_from.get_end())
            arrows_full.append(arrow)

        arrows_partial = []

        for field_line_from, field_line_to in field_lines_partial:
            point = field_line_from.point_from_proportion(2/field_line_from.get_arc_length())
            next_point = field_line_from.point_from_proportion(2/field_line_from.get_arc_length() + 1e-2)
            angle = angle_of_vector(next_point - point)

            arrow_from = ArrowTriangleFilledTip(start_angle=angle, color=WHITE).scale(0.4)
            arrow_from.shift(point - center_of_mass(arrow_from.get_vertices()))

            point = field_line_to.point_from_proportion(1 - 2/field_line_to.get_arc_length())
            next_point = field_line_to.point_from_proportion(1 - 2/field_line_to.get_arc_length() + 1e-2)
            angle = angle_of_vector(next_point - point)

            arrow_to = ArrowTriangleFilledTip(start_angle=angle, color=WHITE).scale(0.4)
            arrow_to.shift(point - center_of_mass(arrow_to.get_vertices()))

            arrows_partial.append((arrow_from, arrow_to))

        self.play(AnimationGroup(
            total_vec_field.animate.set_opacity(0.25),
            AnimationGroup(
                AnimationGroup(*[
                    Succession(
                        Create(field_line_from, rate_func=rate_functions.ease_in_quad),
                        Create(field_line_to, rate_func=rate_functions.ease_out_quad),
                    )
                    for field_line_from, field_line_to in field_lines_full + field_lines_partial
                ]),
                AnimationGroup(
                    *[GrowFromCenter(arrow) for arrow, _ in arrows_partial],
                    rate_func=rate_functions.ease_out_cubic
                ),
                lag_ratio=0.25
            ),
            AnimationGroup(
                AnimationGroup(
                    *[GrowFromCenter(arrow) for arrow in arrows_full],
                    rate_func=rate_functions.ease_out_cubic
                ),
                AnimationGroup(
                    *[GrowFromCenter(arrow) for _, arrow in arrows_partial],
                    rate_func=rate_functions.ease_out_cubic
                ),
                lag_ratio=0.6
            ),
            lag_ratio=0.45, run_time=6
        ))

        self.wait()

        self.play(AnimationGroup(
            FadeOut(total_vec_field),
            run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()
