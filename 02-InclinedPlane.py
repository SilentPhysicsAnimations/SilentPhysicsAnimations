from manim import *

#######################################

GRAVITY_FORCE_LABEL = 'F_g'
GRAVITY_FORCE_COLOR = BLUE

NORMAL_FORCE_LABEL = 'F_N'
NORMAL_FORCE_COLOR = YELLOW

FRICTION_FORCE_LABEL = r'F_{\text{fr}}'
FRICTION_FORCE_COLOR = RED

TOTAL_FORCE_LABEL = r'F_{\text{net}}'
TOTAL_FORCE_COLOR = GREEN

#######################################

TEX = { 'tex_template': TexTemplateLibrary.default }

#######################################

class Block(Rectangle):
    def __init__(self, mass=1, fill_opacity=1, fill_color=WHITE, stroke_width=0, **kwargs):
        self.mass = mass

        Rectangle.__init__(
            self,
            fill_opacity=fill_opacity,
            fill_color=fill_color,
            stroke_width=stroke_width,
            **kwargs
        )

    def get_label(self, label_text, label_constructor=MathTex, height_ratio=0.33, color=BLACK, **kwargs):
        label = label_constructor(label_text, color=color, **kwargs)
        corners = self.get_vertices()
        label.rotate(angle_of_vector(corners[0] - corners[1]))
        label.scale_to_fit_height(height_ratio * self.height)
        label.move_to(self)
        return label

    def get_point_on_block(self, direction):
        corners = self.get_vertices()
        if (direction == ORIGIN).all():
            return sum(corners) / 4
        elif (direction == RIGHT).all():
            return mid(corners[0], corners[3])
        elif (direction == UR).all():
            return corners[0]
        elif (direction == UP).all():
            return mid(corners[0], corners[1])
        elif (direction == UL).all():
            return corners[1]
        elif (direction == LEFT).all():
            return mid(corners[1], corners[2])
        elif (direction == DL).all():
            return corners[2]
        elif (direction == DOWN).all():
            return mid(corners[2], corners[3])
        elif (direction == DR).all():
            return corners[3]
        else:
            raise Exception("Invalid input")


class InclinedPlane(Polygon):
    def __init__(self, angle, max_height=6, max_width=10, color=WHITE, **kwargs):
        self.angle = angle

        height, width = max_height, max_width
        tan = np.tan(angle)
        if width * tan > height:
            width = height / tan
        else:
            height = width * tan

        Polygon.__init__(
            self,
            [0, 0, 0],
            [0, height, 0],
            [width, 0, 0],
            color=color, **kwargs
        )

    def get_angle_arc(self, **kwargs):
        return Arc(
            start_angle=PI-self.angle, angle=self.angle,
            arc_center=self.get_vertices()[2],
            **kwargs
        )

    def get_angle_label(self, label_text, label_contructor=MathTex, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER, radius=1, **kwargs):
        label = label_contructor(label_text, **kwargs)
        label.move_to(
            self.get_vertices()[2] + Circle(radius+buff).point_at_angle(PI - self.angle/2)
        )
        return label


class InclinedPlaneScene(MovingCameraScene):
    def construct(self):
        m = 1
        theta = 30*DEGREES
        mu = 0.25
        g = 9.82

        cos = np.cos(theta)
        sin = np.sin(theta)
        a = g * (sin - mu*cos)

        time = ValueTracker(0)

        plane = InclinedPlane(angle=theta).to_corner(DL)
        block = Block(width=1.618*0.75, height=0.75, mass=m)

        start = plane.get_corner(UL)
        slide_dist = np.sqrt(plane.width**2 + plane.height**2) - block.width 
        slide_time = np.sqrt(2*slide_dist / a)

        block.rotate(-theta, about_point=start)

        def d(t):
            return max(0, 1/2 * g * (sin - mu*cos) * t**2)

        def block_updater(mob):
            t = time.get_value()
            shift_val = mob.get_point_on_block(DL) - mob.get_center()
            mob.move_to(start).shift(d(t) * np.array([cos, -sin, 0]) - shift_val)

        block_updater(block)
        block.add_updater(block_updater)

        self.add(plane, block)

        block_label = block.get_label('m', **TEX)
        self.play(FadeIn(block_label))
        block.add(block_label)

        self.wait(1.75)

        angle_arc = plane.get_angle_arc()
        angle_label = plane.get_angle_label(r'\theta', font_size=35, **TEX)

        self.play(AnimationGroup(
            AnimationGroup(
                Create(angle_arc),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(angle_label, target_position=angle_arc),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.2, run_time=2
        ))

        plane.add(angle_arc, angle_label)

        self.wait(1.25)

        self.play(
            time.animate.set_value(slide_time),
            rate_func=rate_functions.linear,
            run_time=slide_time
        )

        self.wait(2.5)

        half_slide_time = np.sqrt(slide_dist / a)

        self.play(time.animate.set_value(half_slide_time))

        self.wait(2.5)

        forces = {
            'gravity': {
                'size': m*g,
                'angle': -PI/2,
                'origin': DOWN,
                'color': GRAVITY_FORCE_COLOR,
                'label_text': GRAVITY_FORCE_LABEL,
                'label_direction': RIGHT,
            },
            'normal': {
                'size': m*g*cos,
                'angle': PI/2 - theta,
                'origin': UP,
                'color': NORMAL_FORCE_COLOR,
                'label_text': NORMAL_FORCE_LABEL,
                'label_direction': RIGHT,
            },
            'friction': {
                'size': mu*m*g*cos,
                'angle': PI - theta,
                'origin': LEFT,
                'color': FRICTION_FORCE_COLOR,
                'label_text': FRICTION_FORCE_LABEL,
                'label_direction': UL,
            },
            'total': {
                'size': m*a,
                'angle': -theta,
                'origin': RIGHT,
                'color': TOTAL_FORCE_COLOR,
                'label_text': TOTAL_FORCE_LABEL,
                'label_direction': RIGHT,
            }
        }

        vector_config = {
            'tip_length': 0.2,
            'max_tip_length_to_length_ratio': 1,
            'stroke_width': 5,
            'max_stroke_width_to_length_ratio': 1000
        }

        def vector_updater(direction):
            return lambda mob: mob.shift(
                block.get_point_on_block(direction) - mob.get_start()
            )

        def vector_label_updater(vector, direction):
            return lambda mob: mob.next_to(vector.get_tip(), direction)

        each_force = ['gravity', 'normal', 'friction']

        force_scale = 0.2
        for force in [*each_force, 'total']:
            info = forces[force]

            vector = Vector(force_scale*info['size']*RIGHT, color=info['color'], **vector_config)
            vector.rotate(info['angle'], about_point=ORIGIN)
            vector.shift(block.get_point_on_block(info['origin']))

            label = MathTex(info['label_text'], color=info['color'], font_size=35, **TEX)
            label.next_to(vector.get_tip(), info['label_direction'])

            info.update({'vector': vector, 'label': label})
            
            if force == 'total':
                continue

            self.play(AnimationGroup(
                AnimationGroup(
                    GrowArrow(vector),
                    rate_func=rate_functions.ease_out_cubic
                ),
                AnimationGroup(
                    FadeIn(label, target_position=vector.get_tip()),
                    rate_func=rate_functions.ease_out_cubic
                ),
                lag_ratio=0.4, run_time=2
            ))

            self.wait(1.5)

            if force == 'gravity':
                formula = MathTex(info['label_text'], '=', 'mg', font_size=35, **TEX)
                formula.shift(label.get_center() - formula[0].get_center())
                self.play(
                    FadeIn(formula[1:], shift=MED_SMALL_BUFF*RIGHT),
                    rate_func=rate_functions.ease_out_cubic
                )
                self.wait(2)
                self.play(
                    FadeOut(formula[1:], shift=MED_SMALL_BUFF*LEFT),
                    run_time=0.5, rate_func=rate_functions.ease_in_cubic
                )
                self.wait()
            elif force == 'friction':
                formula = MathTex(info['label_text'], '=', r'\mu', forces['normal']['label_text'], font_size=35, **TEX)
                formula.shift(label.get_center() - formula[0].get_center())
                shift_val = formula[-1].get_center() - forces['normal']['label'].get_center()
                self.play(AnimationGroup(
                    AnimationGroup(
                        FadeIn(formula[1:-1], shift=MED_SMALL_BUFF*RIGHT),
                        run_time=1.2, rate_func=rate_functions.ease_out_cubic
                    ),
                    forces['normal']['label'].animate(run_time=3, rate_func=rate_functions.ease_in_out_cubic).shift(shift_val),
                    lag_ratio=0.05
                ))
                self.wait(2)
                self.play(
                    FadeOut(formula[1:-1], shift=MED_SMALL_BUFF*LEFT),
                    forces['normal']['label'].animate.shift(-shift_val),
                    run_time=0.5, rate_func=rate_functions.ease_in_cubic
                )
                self.wait()

        self.play(
            self.camera.frame.animate.shift(2.5*RIGHT),
            run_time=1.25
        )
        
        total_force_eq = MathTex(
            forces['total']['label_text'], 
            '=', forces['gravity']['label_text'], 
            '+', forces['normal']['label_text'], 
            '+', forces['friction']['label_text'],
            **TEX
        )

        total_force_eq.move_to(self.camera.frame, aligned_edge=UR).shift(1.5*DL)
        
        for i, force in zip([0, 2, 4, 6], ['total', *each_force]):
            total_force_eq[i].set(color=forces[force]['color'])

        self.play(AnimationGroup(
            FadeIn(total_force_eq[0:2]),
            AnimationGroup(
                *[
                    forces[force]['label'].animate(run_time=1.5)
                        .move_to(total_force_eq[i])
                        .scale_to_fit_height(total_force_eq[i].height)
                    for force, i in zip(each_force, [2, 4, 6])
                ],
                *[FadeIn(total_force_eq[i]) for i in (3, 5)],
                lag_ratio=0.2
            ),
            lag_ratio=0.2
        ))

        self.remove(*[forces[force]['label'] for force in each_force])
        self.add(total_force_eq[2:7:2])

        self.wait()

        arrow_width = MathTex('F').width
        vec_arrows = {}
        for i, force in zip([0, 2, 4, 6], ['total', *each_force]):
            arrow = Vector(arrow_width*RIGHT, color=forces[force]['color'])
            arrow.next_to(total_force_eq[i], UP, aligned_edge=LEFT, buff=0.1).shift(0.035*RIGHT)
            vec_arrows[force] = arrow

        self.play(AnimationGroup(*[GrowArrow(arr) for arr in vec_arrows.values()], lag_ratio=0.2))

        self.wait(0.2)

        vector_shifts = []
        vector_sum_point = block.get_point_on_block(RIGHT)
        for force in each_force:
            vector = forces[force]['vector']
            vector_shifts.append((vector, vector_sum_point - vector.get_start()))
            vector_sum_point += vector.get_end() - vector.get_start()

        self.play(AnimationGroup(
            *[vector.animate.shift(shift_val) for vector, shift_val in vector_shifts],
            lag_ratio=0.4, run_time=2
        ))
        
        self.wait()

        self.play(AnimationGroup(
            AnimationGroup(
                GrowArrow(forces['total']['vector']),
                total_force_eq[0].animate.move_to(forces['total']['label']).scale_to_fit_height(forces['total']['label'].height),
                lag_ratio=0.2, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeOut(total_force_eq[1:], *vec_arrows.values()),
                FadeOut(*[forces[force]['vector'] for force in each_force]),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.5, run_time=4
        ))

        for vector, shift_val in vector_shifts:
            vector.shift(-shift_val)

        self.remove(total_force_eq[0])
        self.add(forces['total']['label'])

        self.wait(2.5)

        total_vec = forces['total']['vector']

        self.camera.frame.save_state()
        self.play(
            self.camera.frame.animate.move_to(total_vec).set(width=13*total_vec.width),
            run_time=2
        )

        self.wait()

        comp_config = {
            'color': forces['total']['color'],
            'stroke_width': 2,
            'dash_length': 0.03,
        }

        total_x = DashedLine(ORIGIN, total_vec.width*RIGHT, **comp_config).shift(total_vec.get_start())
        total_y = DashedLine(ORIGIN, total_vec.height*DOWN, **comp_config).shift(total_x.get_end())

        self.play(AnimationGroup(
            GrowFromPoint(total_x, total_vec.get_start()),
            GrowFromPoint(total_y, total_x.get_end()),
            lag_ratio=0.5, run_time=2.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(
            FadeOut(forces['total']['label']),
            self.camera.frame.animate.move_to(block).scale(1.5)
        )

        total_x.add_updater(lambda mob: mob.put_start_and_end_on(
            total_vec.get_start(), np.array([total_vec.get_end()[0], total_vec.get_start()[1], 0])
        ))
        total_y.add_updater(lambda mob: mob.put_start_and_end_on(
            np.array([total_vec.get_end()[0], total_vec.get_start()[1], 0]), total_vec.get_end() + 1e-6*DOWN
        ))

        x_vec = Vector(RIGHT, tip_length=0.15, stroke_width=1.25)
        y_vec = x_vec.copy().rotate(PI/2, about_point=ORIGIN)
        x_label = MathTex('x').scale(0.6).next_to(x_vec, RIGHT, buff=SMALL_BUFF)
        y_label = MathTex('y').scale(0.6).next_to(y_vec, UP, buff=SMALL_BUFF)
        axes = VGroup(x_vec, y_vec, x_label, y_label).scale(0.7).move_to(block.get_center()+3*RIGHT)
        self.play(FadeIn(axes))
        self.wait(2)
        self.play(Rotate(axes, -theta, about_point=x_vec.get_start()), run_time=1.5)
        self.wait(2)

        plane.scale(3, about_point=block.get_point_on_block(DOWN))

        block.remove_updater(block_updater)

        rotate_about = block.get_center()

        self.play(
            Rotate(
                Group(plane, block, total_vec, axes),
                theta, about_point=rotate_about
            ),
            run_time=2.5
        )

        for force in each_force:
            forces[force]['vector'].rotate(theta, about_point=rotate_about)

        total_x.clear_updaters()
        total_y.clear_updaters()

        self.play(FadeOut(total_x, total_y), run_time=0.3)

        self.wait(2)

        start = block.get_center()
        time.set_value(0)
        def block_updater_flat(mob):
            t = time.get_value()
            mob.move_to(start + d(t)*RIGHT)

        block.add_updater(block_updater_flat)

        total_vec.add_updater(lambda mob: mob.shift(block.get_point_on_block(RIGHT) - mob.get_start()))

        self.play(
            time.animate.set_value(slide_time),
            run_time=slide_time, rate_func=rate_functions.linear
        )

        self.wait()

        self.play(AnimationGroup(
            self.camera.frame.animate.scale(1.4).shift(3*RIGHT + 0.5*UP),
            FadeOut(axes),
            time.animate.set_value(0),
            lag_ratio=0.2, run_time=5, rate_func=rate_functions.ease_out_cubic
        ))

        total_vec.clear_updaters()

        self.wait(2)

        force_eq = MathTex(forces['total']['label_text'], '=', 'm', 'a', **TEX)
        force_eq.move_to(self.camera.frame, aligned_edge=UR).shift(1.5*DL)
        force_eq[0].set(color=forces['total']['color'])

        vec_arrows['total'].next_to(force_eq[0], UP, aligned_edge=LEFT, buff=0.1).shift(0.035*RIGHT)

        vec_arrows['a'] = Vector(force_eq[3].width*RIGHT, color=WHITE)
        vec_arrows['a'].next_to(force_eq[3], UP, aligned_edge=LEFT, buff=0.1)

        self.play(AnimationGroup(
            ReplacementTransform(total_vec.copy(), vec_arrows['total']),
            *[FadeIn(letter) for letter in force_eq],
            GrowArrow(vec_arrows['a']),
            lag_ratio=0.15, run_time=4, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.wait()

        new_force_eq = MathTex(
            forces['gravity']['label_text'], 
            '+', forces['normal']['label_text'], 
            '+', forces['friction']['label_text'],
            '=', 'm', 'a',
            **TEX
        )

        for i, force in zip([0, 2, 4], each_force):
            new_force_eq[i].set(color=forces[force]['color'])

        new_force_eq.shift(force_eq[-1].get_center() - new_force_eq[-1].get_center())

        shift_val = new_force_eq[0].get_center() - total_force_eq[2].get_center()
        for force in each_force:
            vec_arrows[force].shift(shift_val)

        self.add(new_force_eq[5:])
        force_eq[1:].shift(200*RIGHT)

        self.play(AnimationGroup(
            FadeOut(
                Group(force_eq[0], vec_arrows['total']),
                shift=0.25*UP
            ),
            FadeIn(
                Group(new_force_eq[:5]).add(*[vec_arrows[force] for force in each_force]),
                shift=0.25*UP
            ),
            FadeOut(total_vec),
            *[GrowArrow(forces[force]['vector']) for force in each_force],
            run_time=3, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.wait()

        grav_vec = forces['gravity']['vector']
        grav_vec_corner = np.array([grav_vec.get_start()[0], grav_vec.get_end()[1], 0])
        grav_y_line = DashedLine(grav_vec.get_start(), grav_vec_corner, color=forces['gravity']['color'])
        grav_x_line = DashedLine(grav_vec_corner, grav_vec.get_end(), color=forces['gravity']['color'])
        self.play(AnimationGroup(
            Create(grav_y_line),
            Create(grav_x_line),
            lag_ratio=1, run_time=0.8, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        self.camera.frame.save_state()
        
        to_rotate = Group(
            grav_vec, grav_y_line, grav_x_line,
            forces['normal']['vector'], forces['friction']['vector'],
            block
        )

        not_to_rotate = Group(
            new_force_eq,
            vec_arrows['gravity'], vec_arrows['normal'],
            vec_arrows['friction'], vec_arrows['a']
        )

        shift_val = -self.camera.frame_center

        self.play(
            AnimationGroup(
                plane.animate.rotate(-theta, about_point=rotate_about).scale(1/3, about_point=block.get_point_on_block(DOWN)),
                to_rotate.animate.rotate(-theta, about_point=rotate_about),
                lag_ratio=0.1
            ),
            FadeOut(not_to_rotate, shift=shift_val),
            self.camera.frame.animate.center().set(width=14),
            run_time=3
        )

        self.wait()

        plane.remove(angle_label)
        plane_copy = plane.copy().set(stroke_width=2)
        angle_label_copy = angle_label.copy()
        self.add(angle_label_copy)

        self.play(
            plane_copy.animate
                .flip()
                .shift(grav_vec.get_start() - plane_copy.get_corner(DL))
                .rotate(-PI/2 - theta, about_point=grav_vec.get_start())
                .scale_to_fit_height(grav_vec.get_length(), about_point=grav_vec.get_start()),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        self.play(
            angle_label_copy.animate.scale(0.9).move_to(grav_vec.get_start()).shift(0.3*DOWN),
            AnimationGroup(
                AnimationGroup(
                    to_rotate.animate.rotate(theta, about_point=rotate_about),
                    plane_copy.animate.rotate(theta, about_point=rotate_about).set(stroke_width=0),
                ),
                plane.animate.rotate(theta, about_point=rotate_about).scale(3, about_point=block.get_point_on_block(DOWN)),
                lag_ratio=0.1
            ),
            FadeIn(not_to_rotate, shift=-shift_val),
            Restore(self.camera.frame),
            run_time=3
        )


        self.wait()

        grav_vec_label = MathTex('mg', color=forces['gravity']['color']).scale(0.6).next_to(grav_vec.get_center())

        self.play(
            FadeIn(grav_vec_label, target_position=grav_vec.get_center()),
            rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        grav_x_label = MathTex('mg', r'\sin', r'\theta', color=forces['gravity']['color']).scale(0.6).next_to(grav_x_line, DOWN)
        grav_y_label = MathTex('-', 'mg', r'\cos', r'\theta', color=forces['gravity']['color']).scale(0.6).next_to(grav_y_line, LEFT)

        self.play(
            FadeOut(plane_copy),
            AnimationGroup(
                ReplacementTransform(grav_vec_label, grav_x_label[0]),
                FadeIn(grav_x_label[1]),
                ReplacementTransform(angle_label_copy, grav_x_label[2]),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                ReplacementTransform(grav_vec_label.copy(), grav_y_label[1]),
                FadeIn(grav_y_label[2]),
                ReplacementTransform(angle_label_copy.copy(), grav_y_label[3]),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.3, run_time=3.5
        )

        self.wait()

        self.play(
            Group(new_force_eq[0], vec_arrows['gravity']).animate.shift(1.5*DOWN),
            forces['normal']['vector'].animate.set_opacity(0.25),
            forces['friction']['vector'].animate.set_opacity(0.25)
        )

        equals = MathTex('=').next_to(new_force_eq[0])

        matrix_config = {
            'left_bracket': '(',
            'right_bracket': ')',
            'bracket_v_buff': SMALL_BUFF,
            'v_buff': 0.9,
            'element_to_mobject_config': TEX,
        }

        self.play(FadeIn(grav_y_label[0], shift=0.1*LEFT), rate_func=rate_functions.ease_out_cubic)

        grav_matrix = Matrix(
            [[r'mg\sin\theta'], [r'-mg\cos\theta']],
            **matrix_config
        ).scale(0.7).set(color=forces['gravity']['color'])
        grav_matrix.next_to(equals)

        grav_entries = grav_matrix.get_entries()

        self.play(AnimationGroup(
            FadeIn(Group(equals, grav_matrix.get_brackets()), shift=0.25*RIGHT),
            grav_x_label.animate.scale_to_fit_height(grav_entries[0].height).move_to(grav_entries[0], aligned_edge=RIGHT),
            grav_y_label.animate.scale_to_fit_height(grav_entries[1].height).move_to(grav_entries[1], aligned_edge=RIGHT),
            lag_ratio=0.3, run_time=5, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(grav_x_label, grav_y_label)
        self.add(grav_entries)

        self.wait(1.5)

        self.play(AnimationGroup(
            grav_matrix.animate.move_to(new_force_eq[0], aligned_edge=RIGHT).shift(1.5*UP),
            FadeOut(Group(new_force_eq[0], vec_arrows['gravity'], equals)),
            forces['normal']['vector'].animate.set_opacity(1),
            forces['friction']['vector'].animate.set_opacity(1),
            FadeOut(grav_y_line, grav_x_line),
            run_time=2.5, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.wait()

        self.play(
            Group(new_force_eq[2], vec_arrows['normal']).animate.shift(1.5*DOWN),
            forces['gravity']['vector'].animate.set_opacity(0.25),
            forces['friction']['vector'].animate.set_opacity(0.25)
        )

        norm_matrix = Matrix(
            [['0'], [forces['normal']['label_text']]],
            **matrix_config
        ).scale(0.7).set(color=forces['normal']['color'])

        equals.next_to(new_force_eq[2])
        norm_matrix.next_to(equals)

        self.play(AnimationGroup(
            FadeIn(Group(equals, norm_matrix.get_brackets()), shift=0.25*RIGHT),
            FadeIn(norm_matrix.get_entries()),
            lag_ratio=0.4, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            Group(grav_matrix, new_force_eq[1])
                .animate(run_time=2).shift((norm_matrix.width - new_force_eq[2].width) * LEFT),
            AnimationGroup(
                norm_matrix.animate.move_to(new_force_eq[2], aligned_edge=RIGHT).shift(1.5*UP),
                FadeOut(Group(new_force_eq[2], vec_arrows['normal'], equals)),
                forces['gravity']['vector'].animate.set_opacity(1),
                forces['friction']['vector'].animate.set_opacity(1),
                run_time=2.5, rate_func=rate_functions.ease_in_out_cubic
            ),
            lag_ratio=0.1
        ))

        self.wait()

        fric_matrix = Matrix(
            [[r'-\mu ' + forces['normal']['label_text']], ['0']],
            color=forces['friction']['color'],
            **matrix_config
        ).scale(0.7).set(color=forces['friction']['color'])

        self.play(
            Group(new_force_eq[4], vec_arrows['friction']).animate.shift(1.5*DOWN),
            forces['gravity']['vector'].animate.set_opacity(0.25),
            forces['normal']['vector'].animate.set_opacity(0.25)
        )

        equals.next_to(new_force_eq[4])
        fric_matrix.next_to(equals)

        self.play(AnimationGroup(
            FadeIn(Group(equals, fric_matrix.get_brackets()), shift=0.25*RIGHT),
            FadeIn(fric_matrix.get_entries()),
            lag_ratio=0.4, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            Group(grav_matrix, new_force_eq[1], norm_matrix, new_force_eq[3])
                .animate(run_time=2).shift((fric_matrix.width - new_force_eq[4].width) * LEFT),
            AnimationGroup(
                fric_matrix.animate.move_to(new_force_eq[4], aligned_edge=RIGHT).shift(1.5*UP),
                FadeOut(Group(new_force_eq[4], vec_arrows['friction'], equals)),
                forces['gravity']['vector'].animate.set_opacity(1),
                forces['normal']['vector'].animate.set_opacity(1),
                run_time=2.5, rate_func=rate_functions.ease_in_out_cubic
            ),
            lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        accel_vector = forces['total']['vector'].copy().set(color=WHITE)

        accel_matrix = Matrix(['a', '0'], **matrix_config).scale(0.7)
        accel_matrix.move_to(new_force_eq[7], aligned_edge=LEFT).shift(1.5*DOWN + SMALL_BUFF*RIGHT)
        equals.next_to(accel_matrix, LEFT)

        a_target = new_force_eq[7].copy().next_to(equals, LEFT)

        self.play(
            Group(new_force_eq[7], vec_arrows['a']).animate.shift(a_target.get_center() - new_force_eq[7].get_center()),
            forces['gravity']['vector'].animate.set_opacity(0.25),
            forces['normal']['vector'].animate.set_opacity(0.25),
            forces['friction']['vector'].animate.set_opacity(0.25),
            GrowArrow(accel_vector)
        )

        self.wait()

        self.play(AnimationGroup(
            FadeIn(Group(equals, accel_matrix.get_brackets()), shift=0.25*RIGHT),
            FadeIn(accel_matrix.get_entries()),
            lag_ratio=0.4, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(AnimationGroup(
            accel_matrix.animate.shift(1.5*UP),
            FadeOut(Group(new_force_eq[7], vec_arrows['a'], equals)),
            forces['gravity']['vector'].animate.set_opacity(1),
            forces['normal']['vector'].animate.set_opacity(1),
            forces['friction']['vector'].animate.set_opacity(1),
            FadeOut(accel_vector),
            run_time=2.5, rate_func=rate_functions.ease_in_out_cubic
        ))

        self.wait(2)

        self.play(AnimationGroup(
            FadeOut(block, plane, *[forces[force]['vector'] for force in each_force]),
            self.camera.frame.animate.move_to(Group(grav_matrix, accel_matrix)),
            lag_ratio=0.1, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()
        
        left_plus_up, left_plus_down = new_force_eq[1], new_force_eq[1].copy()
        right_plus_up, right_plus_down = new_force_eq[3], new_force_eq[3].copy()
        equals_up, equals_down = new_force_eq[5], new_force_eq[5].copy()
        m_sign = new_force_eq[6]

        y_up = fric_matrix.get_entries()[0].get_y()
        y_down = norm_matrix.get_entries()[1].get_y()

        brace = Brace(grav_matrix, LEFT)

        self.play(AnimationGroup(
            AnimationGroup(
                GrowFromCenter(brace),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeOut(Group(*[matrix.get_brackets() for matrix in [grav_matrix, norm_matrix, fric_matrix, accel_matrix]])),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                left_plus_up.animate.scale(0.7).set_y(y_up),
                left_plus_down.animate.scale(0.7).set_y(y_down),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                right_plus_up.animate.scale(0.7).set_y(y_up),
                right_plus_down.animate.scale(0.7).set_y(y_down),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                equals_up.animate.scale(0.7).set_y(y_up),
                equals_down.animate.scale(0.7).set_y(y_down),
                rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                m_sign.animate().scale(0.7).next_to(accel_matrix.get_entries()[0], LEFT, buff=0.03),
                rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.05, run_time=4, 
        ))

        self.wait(0.6)

        x_eq = MathTex(
            'm', 'g', r'\sin', r'\theta', '-', r'\mu', forces['normal']['label_text'], '=', 'm', 'a',
            **TEX
        ).scale(0.7)
        y_eq = MathTex(
            '-', 'm', 'g', r'\cos', r'\theta', '+', forces['normal']['label_text'], '=', '0',
            **TEX
        ).scale(0.7)

        Group(x_eq, y_eq).arrange_in_grid(2, 1, col_alignments='l').move_to(self.camera.frame_center)
        x_eq.set_y(y_up)
        y_eq.set_y(y_down)

        x_eq[:4].set(color=forces['gravity']['color'])
        y_eq[1:5].set(color=forces['gravity']['color'])
        y_eq[6].set(color=forces['normal']['color'])
        x_eq[5:7].set(color=forces['friction']['color'])

        grav_matrix_y = MathTex('-', r'mg\cos\theta').set(color=forces['gravity']['color'])
        grav_matrix_y.scale(0.7).move_to(grav_matrix.get_entries()[1])
        fric_matrix_x = MathTex(
            '-', r'\mu ' + forces['normal']['label_text'],
            **TEX
        ).set(color=forces['friction']['color'])
        fric_matrix_x.scale(0.7).move_to(fric_matrix.get_entries()[0])

        self.remove(grav_matrix.get_entries()[1], fric_matrix.get_entries()[0])
        self.add(grav_matrix_y, fric_matrix_x)

        replacement_map = {
            grav_matrix.get_entries()[0]: x_eq[:4],
            grav_matrix_y[1]: y_eq[1:5],
            norm_matrix.get_entries()[1]: y_eq[6],
            fric_matrix_x[1]: x_eq[5:7],
            accel_matrix.get_entries()[0]: x_eq[9],
            accel_matrix.get_entries()[1]: y_eq[8],
            left_plus_down: y_eq[5],
            equals_up: x_eq[7],
            equals_down: y_eq[7],
            m_sign: x_eq[8],
        }

        self.play(AnimationGroup(
            AnimationGroup(
                FadeOut(left_plus_up),
                FadeOut(Group(right_plus_up, right_plus_down)),
                FadeOut(norm_matrix.get_entries()[0]),
                FadeOut(fric_matrix.get_entries()[1]),
                lag_ratio=0.1
            ),
            AnimationGroup(
                *[old.animate.move_to(new) for old, new in replacement_map.items()],
                grav_matrix_y[0].animate.move_to(y_eq[0]).set(color=WHITE),
                fric_matrix_x[0].animate.move_to(x_eq[4]).set(color=WHITE),
                brace.animate.next_to(Group(x_eq, y_eq), LEFT),
                run_time=1.5
            ),
            lag_ratio=0.7,
        ))

        self.remove(*replacement_map)
        self.remove(grav_matrix_y[0], fric_matrix_x[0])
        self.add(x_eq, y_eq)

        self.wait()

        self.camera.frame.save_state()
        camera_width = 3 * y_eq.width
        self.play(AnimationGroup(
            self.camera.frame.animate.set_width(camera_width).move_to(y_eq),
            x_eq.animate.set_opacity(0.25),
            brace.animate.set_opacity(0.25),
            run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        flipped_y_eq = MathTex(
            forces['normal']['label_text'], '=', 'm', 'g', r'\cos', r'\theta',
            **TEX
        ).scale(0.7)
        flipped_y_eq.shift(y_eq[6].get_center() - flipped_y_eq[0].get_center())

        self.wait(2)

        self.play(AnimationGroup(
            *[FadeOut(y_eq[i]) for i in [0, 5, 8]],
            MoveAlongPath(
                y_eq[1:5],
                ArcBetweenPoints(
                    y_eq[1:5].get_center(),
                    flipped_y_eq[2:].get_center()
                )
            ),
            lag_ratio=0.2, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        shift_val = (y_eq[6].get_left()[0] - x_eq.get_left()[0]) * LEFT
        self.play(
            y_eq[6:8].animate.shift(shift_val),
            y_eq[1:5].animate.shift(shift_val).set(color=forces['normal']['color']),
            Restore(self.camera.frame),
            x_eq.animate.set_opacity(1),
            brace.animate.set_opacity(1),
            run_time=3, rate_func=rate_functions.ease_in_out_cubic
        )

        self.wait()

        combined_eq = MathTex(
            'm', 'g', r'\sin', r'\theta', '-', r'\mu', 'm', 'g', r'\cos', r'\theta', '=', 'm', 'a',
            **TEX
        ).scale(0.7)
        combined_eq.shift(x_eq[0].get_center() - combined_eq[0].get_center())

        self.play(AnimationGroup(
            AnimationGroup(
                *[x_eq[i].animate.move_to(combined_eq[i+3]) for i in [7, 8, 9]],
            ),
            FadeOut(x_eq[6], shift=0.5*UP),
            y_eq[1:5].animate.move_to(combined_eq[6:10]),
            AnimationGroup(
                self.camera.frame.animate.set_width(camera_width).move_to(combined_eq),
                FadeOut(Group(brace, y_eq[6:8]))
            ),
            lag_ratio=0.05, run_time=5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        self.play(
            FadeOut(Group(x_eq[0], x_eq[8], y_eq[1]), shift=0.5*DOWN),
            run_time=1.5, rate_func=rate_functions.ease_in_out_cubic
        )

        self.wait(0.5)

        final_eq = MathTex(
            'a', '=', 'g', r'\sin', r'\theta', '-', r'\mu', 'g', r'\cos', r'\theta',
            **TEX
        ).scale(0.7).move_to(self.camera.frame_center)

        self.play(
            x_eq[1:6].animate.move_to(final_eq[2:7]),
            y_eq[2:5].animate.move_to(final_eq[7:]),
            MoveAlongPath(x_eq[9], ArcBetweenPoints(x_eq[9].get_center(), final_eq[0].get_center())),
            MoveAlongPath(x_eq[7], ArcBetweenPoints(x_eq[7].get_center(), final_eq[1].get_center())),
            run_time=3, rate_func=rate_functions.ease_in_out_cubic
        )

        self.wait()

