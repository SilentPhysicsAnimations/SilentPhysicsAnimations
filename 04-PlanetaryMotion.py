from manim import *

########################

STAR_COLOR = YELLOW
MAIN_PLANET_COLOR = BLUE

SMALL_PLANET_COLOR = RED
BIG_PLANET_COLOR = GREEN

########################

class Sweep(VMobject):
    def __init__(self, center, func, **kwargs):
        super().__init__(**kwargs)
        self.start_new_path(center)
        self.add_line_to(func())
        self.add_line_to(center)
        self.add_updater(lambda mob:
            mob.set_points(mob.points[:-4]).add_line_to(func()).add_line_to(center)
        )

class PlanetaryMotionScene(MovingCameraScene):
    def construct(self):
        G, m, M = 1, 1, 80
        a, b = 4, 3
        c = np.sqrt(a*a - b*b)

        T = 2 * PI * np.sqrt(a**3 / (G*M))

        star = Dot(radius=0.3, color=STAR_COLOR).shift(c * LEFT)
        planet = Dot(radius=0.1, color=MAIN_PLANET_COLOR).shift(a * RIGHT)

        time_scale = ValueTracker(1)

        def planet_updater(v_0):
            v = v_0

            def updater(mob, dt):
                nonlocal v
                dt *= time_scale.get_value()
                if dt <= 0:
                    return
                n = 10
                dt /= n
                for _ in range(n):
                    r = mob.get_center() - star.get_center()
                    r_sz = np.linalg.norm(r)
                    gravity = G * m * M / (r_sz**2)
                    v -= gravity * dt * r / r_sz
                    mob.shift(v * dt)

            return updater

        planet.add_updater(planet_updater(np.sqrt(G*M * (2/(a + c) - 1/a)) * UP))

        orbit_trace = TracedPath(planet.get_center, stroke_color=GRAY, stroke_width=2, sheen_factor=-0.25, sheen_direction=RIGHT)

        self.add(orbit_trace, star, planet)

        self.wait(T)

        orbit = Ellipse(width=2*a, height=2*b, color=GRAY, stroke_width=2, sheen_factor=-0.25, sheen_direction=RIGHT)

        self.remove(orbit_trace)
        self.add(orbit, star, planet)

        self.wait(3*T/4)
        self.play(
            time_scale.animate.set_value(0),
            run_time=T/2, rate_func=rate_functions.linear
        )

        other_planet = planet.copy()
        other_planet.clear_updaters()

        self.wait(1.5)

        self.play(
            MoveAlongPath(other_planet, orbit.get_subcurve(0, 0.5)),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        other_planet.add_updater(planet_updater(np.sqrt(G*M * (2/(a - c) - 1/a)) * DOWN))

        sweep = Sweep(
            star.get_center(), planet.get_center,
            color=MAIN_PLANET_COLOR, stroke_width=2.5,
            fill_opacity=0.15, stroke_opacity=0.5,
        )

        other_sweep = Sweep(
            star.get_center(), other_planet.get_center,
            color=MAIN_PLANET_COLOR, stroke_width=2.5,
            fill_opacity=0.15, stroke_opacity=0.5,
        )

        self.add(sweep, other_sweep, star, planet)
        time_scale.set_value(0.5)

        self.wait(T/5)

        sweep.clear_updaters()
        other_sweep.clear_updaters()

        self.wait(T/5)

        A_1_label = MathTex("A_1", color=MAIN_PLANET_COLOR).scale(0.8)
        A_1_label.move_to(center_of_mass([other_sweep.points[0], other_sweep.points[-4], other_sweep.points[4], other_sweep.points[len(other_sweep.points) // 2]]))
        A_2_label = MathTex("A_2", color=MAIN_PLANET_COLOR).scale(0.8)
        A_2_label.move_to(center_of_mass([sweep.points[0], sweep.points[-4], sweep.points[4], sweep.points[len(sweep.points) // 2]]))

        self.play(AnimationGroup(
            FadeIn(A_1_label, shift=0.2*UP),
            FadeIn(A_2_label, shift=0.2*UP),
            run_time=2*T/5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(T/5)

        A_eq = MathTex('A_1', '=', 'A_2', color=MAIN_PLANET_COLOR).scale(0.8).next_to(orbit, DOWN, buff=0.5)
        A_eq[1].set(color=WHITE)

        self.play(AnimationGroup(
            ReplacementTransform(A_1_label, A_eq[0]),
            ReplacementTransform(A_2_label, A_eq[2]),
            FadeIn(A_eq[1], shift=0.1*UP),
            self.camera.frame.animate.move_to(Group(orbit, A_eq)),
            run_time=3*T/4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(3*T/4)

        other_planet_opacity = ValueTracker(1)
        other_planet.add_updater(lambda mob: mob.set_opacity(other_planet_opacity.get_value()))

        self.play(AnimationGroup(
            FadeOut(A_eq, shift=0.1*DOWN),
            FadeOut(sweep),
            FadeOut(other_sweep),
            star.animate.set_opacity(1), # Keep star in foreground
            other_planet_opacity.animate.set_value(0),
            self.camera.frame.animate.center(),
            run_time=3*T/4, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(other_planet)

        self.wait(T/2)

        self.play(time_scale.animate.set_value(0), run_time=T/2, rate_func=rate_functions.linear)
        planet.clear_updaters()
        
        self.play(AnimationGroup(
            star.animate.center(),
            orbit.animate.stretch_to_fit_width(2*b).set_sheen(-0.05),
            planet.animate.move_to(b*RIGHT),
            run_time=4.5, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        r_line = DashedLine(star.get_right(), planet.get_center(), color=MAIN_PLANET_COLOR)
        r_label = MathTex('r', color=MAIN_PLANET_COLOR).scale(0.8).next_to(r_line, DOWN, buff=SMALL_BUFF)

        M_label = MathTex('M', color=STAR_COLOR).scale(0.8).next_to(star, LEFT, buff=SMALL_BUFF)
        m_label = MathTex('m', color=MAIN_PLANET_COLOR).scale(0.8).next_to(planet, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            AnimationGroup(
                Create(r_line),
                FadeIn(r_label, shift=SMALL_BUFF*DOWN),
                lag_ratio=0.2, rate_func=rate_functions.ease_out_cubic
            ),
            AnimationGroup(
                FadeIn(M_label, shift=SMALL_BUFF*LEFT),
                FadeIn(m_label, shift=SMALL_BUFF*RIGHT),
                lag_ratio=0.1, rate_func=rate_functions.ease_out_cubic
            ),
            lag_ratio=0.5, run_time=5
        ))

        self.wait(2)

        grav_eq = MathTex('F_G', '=', 'G', '{m', 'M', r'\over', 'r', '^2}').scale(0.8)
        grav_eq.next_to(orbit, buff=LARGE_BUFF).align_to(orbit, UP)
        grav_eq[0:5:4].set(color=STAR_COLOR)
        grav_eq[3:7:3].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            self.camera.frame.animate.move_to(Group(orbit, grav_eq)),
            FadeIn(grav_eq[0:3], shift=SMALL_BUFF*RIGHT),
            ReplacementTransform(m_label.copy(), grav_eq[3]),
            ReplacementTransform(M_label.copy(), grav_eq[4]),
            FadeIn(grav_eq[5], shift=SMALL_BUFF*RIGHT),
            ReplacementTransform(r_label.copy(), grav_eq[6]),
            FadeIn(grav_eq[7], shift=SMALL_BUFF*RIGHT),
            lag_ratio=0.05, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        v_vec = Vector(UP, color=MAIN_PLANET_COLOR).shift(planet.get_center())
        v_label = MathTex('v', color=MAIN_PLANET_COLOR).scale(0.8).next_to(v_vec, buff=SMALL_BUFF)

        self.play(AnimationGroup(
            GrowArrow(v_vec),
            FadeIn(v_label, shift=SMALL_BUFF*RIGHT),
            lag_ratio=0.2, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        centripedal_eq = MathTex('F_c', '=', 'm', '{v', '^2', r'\over', 'r}').scale(0.8)
        centripedal_eq.next_to(grav_eq, DOWN, buff=MED_LARGE_BUFF).align_to(grav_eq, LEFT)
        centripedal_eq[0::3].set(color=MAIN_PLANET_COLOR)
        centripedal_eq[2].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            FadeIn(centripedal_eq[0:2], shift=SMALL_BUFF*RIGHT),
            ReplacementTransform(m_label.copy(), centripedal_eq[2]),
            ReplacementTransform(v_label.copy(), centripedal_eq[3]),
            FadeIn(centripedal_eq[4:6], shift=SMALL_BUFF*RIGHT),
            ReplacementTransform(r_label.copy(), centripedal_eq[6]),
            lag_ratio=0.05, run_time=6, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        planet.add_updater(planet_updater(np.sqrt(G * M / b) * UP))

        self.play(AnimationGroup(
            FadeOut(Group(r_line, r_label, M_label, m_label, v_label, v_vec), rate_func=rate_functions.ease_out_cubic),
            time_scale.animate(rate_func=rate_functions.linear).set_value(1),
            lag_ratio=0.2, run_time=2
        ))

        self.wait(4)

        force_eq = MathTex('F_G', '=', 'F_c').scale(0.8).next_to(centripedal_eq, DOWN, buff=LARGE_BUFF)
        force_eq[0].set(color=STAR_COLOR)
        force_eq[2].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(grav_eq[0].copy(), force_eq[0]),
            FadeIn(force_eq[1], shift=SMALL_BUFF*UP),
            ReplacementTransform(centripedal_eq[0].copy(), force_eq[2]),
            lag_ratio=0.05, run_time=2, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(3)

        force_eq_2 = MathTex('G', '{m', 'M', r'\over', 'r', '^2}', '=', 'm', '{v', '^2', r'\over', 'r}').scale(0.8)
        force_eq_2.next_to(orbit, buff=LARGE_BUFF)
        force_eq_2[1:8:3].set(color=MAIN_PLANET_COLOR)
        force_eq_2[8::3].set(color=MAIN_PLANET_COLOR)
        force_eq_2[2].set(color=STAR_COLOR)

        self.play(AnimationGroup(
            FadeOut(force_eq[0]),
            FadeOut(grav_eq[:2]),
            ReplacementTransform(grav_eq[2:], force_eq_2[:6]),
            ReplacementTransform(force_eq[1], force_eq_2[6]),
            FadeOut(force_eq[2]),
            FadeOut(centripedal_eq[:2]),
            ReplacementTransform(centripedal_eq[2:], force_eq_2[7:]),
            lag_ratio=0.02, run_time=4, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        force_eq_3 = MathTex('G', '{M', r'\over', 'r}', '=', 'v', '^2').scale(0.8)
        force_eq_3.next_to(orbit, buff=LARGE_BUFF)
        force_eq_3[3:6:2].set(color=MAIN_PLANET_COLOR)
        force_eq_3[1].set(color=STAR_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(force_eq_2[0], force_eq_3[0]),
            FadeOut(force_eq_2[1], shift=SMALL_BUFF*LEFT),
            ReplacementTransform(force_eq_2[2], force_eq_3[1]),
            ReplacementTransform(force_eq_2[3], force_eq_3[2]),
            ReplacementTransform(force_eq_2[4], force_eq_3[3]),
            FadeOut(force_eq_2[5], shift=SMALL_BUFF*DOWN),
            ReplacementTransform(force_eq_2[6], force_eq_3[4]),
            FadeOut(force_eq_2[7], shift=SMALL_BUFF*LEFT),
            ReplacementTransform(force_eq_2[8], force_eq_3[5]),
            FadeOut(force_eq_2[11], shift=SMALL_BUFF*DOWN),
            FadeOut(force_eq_2[10], shift=SMALL_BUFF*DOWN),
            ReplacementTransform(force_eq_2[9], force_eq_3[6]),
            lag_ratio=0.02, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        v_eq = MathTex('v', '=', '{2', r'\pi', 'r', r'\over', 'T}').scale(0.8)
        v_eq.next_to(force_eq_3, DOWN, buff=MED_LARGE_BUFF)
        v_eq[0].set(color=MAIN_PLANET_COLOR)
        v_eq[4::2].set(color=MAIN_PLANET_COLOR)

        self.play(
            FadeIn(v_eq, shift=MED_SMALL_BUFF*RIGHT),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait()

        force_eq_4 = MathTex('G', '{M', r'\over', 'r}', '=', r'\left(', '{2', r'\pi', 'r', r'\over', 'T}', r'\right)', '^2').scale(0.8)
        force_eq_4.next_to(orbit, buff=LARGE_BUFF)
        force_eq_4[1].set(color=STAR_COLOR)
        force_eq_4[3].set(color=MAIN_PLANET_COLOR)
        force_eq_4[8:11:2].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(force_eq_3[:5], force_eq_4[:5]),
            ReplacementTransform(force_eq_3[6], force_eq_4[12]),
            FadeOut(force_eq_3[5], shift=SMALL_BUFF*UP),
            FadeIn(force_eq_4[5], shift=SMALL_BUFF*UP),
            FadeIn(force_eq_4[11], shift=SMALL_BUFF*UP),
            ReplacementTransform(v_eq[2:], force_eq_4[6:11]),
            FadeOut(v_eq[:2], shift=SMALL_BUFF*DOWN),
            lag_ratio=0.02, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait()

        force_eq_5 = MathTex('G', '{M', r'\over', 'r}', '=', '{4', r'\pi', '^2', 'r', '^2', r'\over', 'T', '^2}').scale(0.8)
        force_eq_5.next_to(orbit, buff=LARGE_BUFF)
        force_eq_5[1].set(color=STAR_COLOR)
        force_eq_5[3].set(color=MAIN_PLANET_COLOR)
        force_eq_5[8:12:3].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(force_eq_4[:5], force_eq_5[:5]),
            FadeOut(force_eq_4[5], shift=SMALL_BUFF*LEFT),
            FadeOut(force_eq_4[11], shift=SMALL_BUFF*LEFT),
            ReplacementTransform(force_eq_4[9], force_eq_5[10]),
            ReplacementTransform(force_eq_4[6], force_eq_5[5]),
            ReplacementTransform(force_eq_4[7], force_eq_5[6]),
            ReplacementTransform(force_eq_4[8], force_eq_5[8]),
            ReplacementTransform(force_eq_4[10], force_eq_5[11]),
            ReplacementTransform(force_eq_4[12].copy(), force_eq_5[7]),
            ReplacementTransform(force_eq_4[12].copy(), force_eq_5[9]),
            ReplacementTransform(force_eq_4[12], force_eq_5[12]),
            lag_ratio=0.02, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.wait(2)

        force_eq_6 = MathTex('{G', 'M', r'\over', '4', r'\pi', '^2}', '=', '{r', '^3', r'\over', 'T', '^2}').scale(0.8)
        force_eq_6.next_to(orbit, buff=LARGE_BUFF)
        force_eq_6[1].set(color=STAR_COLOR)
        force_eq_6[7::3].set(color=MAIN_PLANET_COLOR)

        self.play(AnimationGroup(
            ReplacementTransform(force_eq_5[0], force_eq_6[0]),
            ReplacementTransform(force_eq_5[1], force_eq_6[1]),
            ReplacementTransform(force_eq_5[2], force_eq_6[2]),
            Transform(force_eq_5[3], force_eq_6[7]),
            ReplacementTransform(force_eq_5[4], force_eq_6[6]),
            ReplacementTransform(force_eq_5[5:8], force_eq_6[3:6]),
            ReplacementTransform(force_eq_5[9], force_eq_6[8]),
            ReplacementTransform(force_eq_5[8], force_eq_6[7]),
            ReplacementTransform(force_eq_5[10], force_eq_6[9]),
            ReplacementTransform(force_eq_5[11:], force_eq_6[10:]),
            lag_ratio=0.02, run_time=3, rate_func=rate_functions.ease_out_cubic
        ))

        self.remove(force_eq_5[3])

        self.wait(5)

        r_2 = 1
        orbit_2 = orbit.copy().scale(r_2/b)
        planet_2 = planet.copy().set(color=SMALL_PLANET_COLOR).scale(2/3).clear_updaters().move_to(orbit_2.get_right())
        planet_2_opacity = ValueTracker(0)
        planet_2.add_updater(lambda mob: mob.set_opacity(planet_2_opacity.get_value()), call_updater=True)
        orbit_2.add_updater(lambda mob: mob.set_stroke(opacity=planet_2_opacity.get_value()), call_updater=True)
        planet_2.add_updater(planet_updater(np.sqrt(G * M / r_2) * UP))

        self.add(orbit_2, planet_2)

        self.play(
            planet_2_opacity.animate.set_value(1),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        planet_2_eq = MathTex('=', '{r', '^3', r'\over', 'T', '^2}').scale(0.8)
        planet_2_eq.next_to(force_eq_6, buff=planet_2_eq[1].get_left()[0]-planet_2_eq[0].get_right()[0])
        planet_2_eq[1:5:3].set(color=SMALL_PLANET_COLOR)

        self.play(
            FadeIn(planet_2_eq, shift=MED_SMALL_BUFF*LEFT),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(2)

        r_3 = 2
        orbit_3 = orbit.copy().scale(r_3/b)
        planet_3 = planet.copy().set(color=BIG_PLANET_COLOR).scale(3/2).clear_updaters().move_to(orbit_3.get_right())
        planet_3_opacity = ValueTracker(0)
        planet_3.add_updater(lambda mob: mob.set_opacity(planet_3_opacity.get_value()), call_updater=True)
        orbit_3.add_updater(lambda mob: mob.set_stroke(opacity=planet_3_opacity.get_value()), call_updater=True)
        planet_3.add_updater(planet_updater(np.sqrt(G * M / r_3) * UP))

        self.add(orbit_3, planet_3)
        
        self.play(
            planet_3_opacity.animate.set_value(1),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        planet_3_eq = MathTex('=', '{r', '^3', r'\over', 'T', '^2}').scale(0.8)
        planet_3_eq.next_to(planet_2_eq, buff=planet_3_eq[1].get_left()[0]-planet_3_eq[0].get_right()[0])
        planet_3_eq[1:5:3].set(color=BIG_PLANET_COLOR)

        self.play(
            FadeIn(planet_3_eq, shift=MED_SMALL_BUFF*LEFT),
            run_time=3, rate_func=rate_functions.ease_out_cubic
        )

        self.wait(4)
