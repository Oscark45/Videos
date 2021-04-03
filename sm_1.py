# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 20:29:33 2021

@author: Oscar
"""

from manim import *
from manimlib.imports import *

from manimlib.utils.rate_functions import linear, slow_into, inverse_square


class SM_intro1(Scene):
    
    # 1. Creating the intro
    
    def construct(self):
        
        course_name = TextMobject(
            "PHYS-A0130 Sähkömagnetismi",
            tex_to_color_map={"text": YELLOW}
        )
        
        # Redundant
        example_text = TexMobject(
            "\\sum_{k=1}^\\infty {1 \\over k^2} = {\\pi^2 \\over 6}",
        )
        vko_aihe = TextMobject(
            "Viikko 1, Sähköstatiikka 1/2",
        )
        group = VGroup(course_name) # Keeping here for reference
        group.arrange(DOWN)
        group.set_width(FRAME_WIDTH - 2 * LARGE_BUFF)
        group.move_to(np.array([0,2.5,0]))

        self.play(FadeIn(course_name))
        self.play(Write(vko_aihe))
        self.wait()


class ChargeField(Scene):
    
    # 2. Point charge electric field.
    
    def construct(self):
        func = lambda p: np.array([
            p[0]/2,  # x
            p[1]/2,  # y
            0        # z
        ])
        
        dot1 = Dot([0,0,0], color=RED, radius=0.35)
        Q_text = TexMobject('+q', color=WHITE, point = [20,0,0])
        self.add(dot1, Q_text)
        
        q_field = VectorField(func, length_func=inverse_square, delta_x=1,delta_y=1)
        self.wait(1)
        self.play(*[GrowArrow(vec) for vec in q_field])
        self.wait(2)
        
        Q_text2 = TexMobject('-q', color=WHITE)
        q_field2 = VectorField(func, length_func=inverse_square2, delta_x=1,delta_y=1)
        
        self.play(ReplacementTransform(Q_text, Q_text2), \
                  ReplacementTransform(q_field, q_field2))

        self.wait(2)
        E_text = TexMobject("\\boldsymbol{E} = {q \\over 4\\pi r^2} \\hat{r}",)
        E_text.to_edge(LEFT)
        self.play(FadeIn(E_text))
        
class CoordDot(GraphScene):
    
    # 3. Calculating the electric field of a charge distribution.
    
    def construct(self):
        
        # Construct origin
        origin = Dot([-4,-3,0], color=WHITE, radius=0.1)
        O_text = TexMobject("0", point = [0,0,0])
        O_text.move_to(point_or_mobject = [-4.5,-3,0])
        self.add(origin, O_text)
        
        # Set electric field position
        E_pos = Arrow((-4,-3,0), (2,0,0), stroke_width = 2)
        E_text = TexMobject("\\boldsymbol{E}(\\boldsymbol{r}) = ?",)
        #E_text.to_edge(LEFT)
        E_text.move_to(point_or_mobject = (2.8,0.5,0))
        r_vec = TexMobject("\\boldsymbol{r}")
        r_vec.move_to(point_or_mobject = (-1, -2, 0))
        
        self.play(ShowCreation(E_pos), FadeIn(r_vec))
        self.play(FadeIn(E_text))
        
        # Set few discrete charges
        q1 = Arrow((-4,-3,0), (-5,2,0), stroke_width = 1)
        q1_dot = Dot((-5,2,0), color=RED, radius=0.1)
        q1_text = TexMobject("q_1")
        q1_text.move_to(point_or_mobject = (-5,2.4,0))
        q1_vec = TexMobject("r_1")
        q1_vec.move_to(point_or_mobject = (-5,-2,0))
        
        q2 = Arrow((-4,-3,0), (-3,1,0), stroke_width = 1)
        q2_dot = Dot((-3,1,0), color=RED, radius=0.1)
        q2_text = TexMobject("q_2")
        q2_text.move_to(point_or_mobject = (-3,1.3,0))
        q2_vec = TexMobject("r_2")
        q2_vec.move_to(point_or_mobject = (-3.8,-0.8,0))
        
        q3 = Arrow((-4,-3,0), (-2,0,0), stroke_width = 1)
        q3_dot = Dot((-2,0,0), color=RED, radius=0.1)
        q3_text = TexMobject("q_3")
        q3_text.move_to(point_or_mobject = (-2,0.3,0))
        q3_vec = TexMobject("r_3")
        q3_vec.move_to(point_or_mobject = (-2.5,-1.5,0))
        
        title_disc = TextMobject("Diskreetti")
        title_disc.move_to(np.array([0,3.2,0]))
        title_cont = TextMobject("Jatkuva")
        title_cont.move_to(np.array([0,3.2,0]))
        
        self.play(ShowCreation(q1), ShowCreation(q2), ShowCreation(q3), \
                  FadeIn(q1_vec), FadeIn(q2_vec), FadeIn(q3_vec),
                  FadeIn(q1_text), FadeIn(q2_text), FadeIn(q3_text), \
                  FadeIn(q1_dot), FadeIn(q2_dot), FadeIn(q3_dot), FadeIn(title_disc),
                  )
        
        E_disc = TexMobject("\\boldsymbol{E}(\\boldsymbol{r})"+\
                            " = \\dfrac{1}{4\\pi \\epsilon_0}"+\
                                "\\sum_i \\dfrac{\\boldsymbol{r}-\\boldsymbol{r}_i'}"+\
                                    "{|\\boldsymbol{r}-\\boldsymbol{r}_i'|^3} q_i ",)
        E_disc.move_to(E_text)
        
        # Continuous distribution
        
        E_cont = TexMobject("E(\\boldsymbol{r}) = ? ",)
        E_cont.move_to(E_text)
        #E_cont.set_y(-1)
        
        Q = Dot((-3,1,0), color=RED, radius=0.8)
        Q_vec = TexMobject("r_i'")
        Q_vec.move_to(point_or_mobject=q2_vec)
        Q_arrow = Arrow((-4,-3,0), (-3,0.4,0), stroke_width = 1)
        
        Q_sym = TexMobject("Q'")
        Q_sym.move_to(point_or_mobject=q2_text)
        Q_sym.set_y(2.3)
        
        #self.play(ReplacementTransform(E_text, E_disc))
        self.play(FadeOut(E_text))
        self.play(FadeIn(E_disc))
        self.play(FadeOut(q1), FadeOut(q3), \
                  FadeOut(q1_text), FadeOut(q3_text), \
                FadeOut(q1_vec), FadeOut(q3_vec),\
                    FadeOut(q1_dot), FadeOut(q3_dot),
                    ReplacementTransform(q2_dot, Q),
                    ReplacementTransform(q2, Q_arrow),
                    ReplacementTransform(q2_vec, Q_vec),
                    ReplacementTransform(q2_text, Q_sym), 
                    FadeOut(E_disc),
                    FadeOut(title_disc),
                    FadeIn(title_cont)
                    )
        self.play(FadeIn(E_cont))
        E_cont2 = TexMobject("E(\\boldsymbol{r}) = \\dfrac{1}{4\\pi \\epsilon_0} \\int  " + \
                             "\\dfrac{\\boldsymbol{r}-\\boldsymbol{r}_i'}" +\
                                    "{|\\boldsymbol{r}-\\boldsymbol{r}_i'|^3} dq'")
        E_cont2.move_to(E_cont)
        self.play(FadeOut(E_cont), FadeIn(E_cont2))
        rs_vec = Arrow((-3,0,0), (1.5,-0.2,0), stroke_width = 1)
        rs_text = TexMobject("\\boldsymbol{r}_s")
        rs_text.move_to((-0.8, -0.4, 0))
        
        E_cont_rs = TexMobject("= \\dfrac{1}{4\\pi \\epsilon_0} \\int  " + \
                             "\\dfrac{d q'}" +\
                                    "{\\boldsymbol{r}_s^2} \\hat{\\boldsymbol{r}}_s}")
        E_cont_rs.next_to(E_cont2, np.array((0., -1.5, 0.)))
        self.play(FadeIn(E_cont_rs), ShowCreation(rs_vec), FadeIn(rs_text))

        E_approx = TexMobject("\\approx \\dfrac{1}{4\\pi \\epsilon_0} \\sum_i  " + \
                             "\\dfrac{\\Delta q'}" +\
                                    "{\\boldsymbol{r}_s^2} \\hat{\\boldsymbol{r}}_s}")
        E_approx.next_to(E_cont_rs, np.array((0., -1.5, 0.)))
        self.play(FadeIn(E_approx) )


class rho_confs(Scene):
    
    # 4. 1D, 2D and 3D charge distributions.
    
    def construct(self):
        
        dim1_text = TextMobject("1D:")
        dim2_text = TextMobject("2D:")
        dim3_text = TextMobject("3D:")
        
        dim1_text.move_to((-6,3,0))
        dim2_text.move_to((-3.5,3,0))
        dim3_text.move_to((-1,3,0))
        
        color = BLUE
        opacity = 0.3
        # 1D
        dom1 = Rectangle(height = 3, width = 0.1, fill_opacity = opacity, color=color)
        dom1.next_to(dim1_text, np.array([0,-3,0]))
        dl_elem = Rectangle(height = 0.2, width = 0.1, fill_opacity = 1, color='white')
        dl_elem.move_to(dom1)
        E1 = TexMobject("E(\\boldsymbol{r}) = \\dfrac{1}{4\\pi \\epsilon_0} \\int_{l'}  " + \
                             "\\dfrac{\\lambda dl}" +\
                                    "{\\boldsymbol{r}_s^2} \\hat{\\boldsymbol{r}}_s}")    
        
            
        lamb = TexMobject("dl")
        lamb.next_to(dl_elem)
        dl = TexMobject("dq = \\lambda dl")
        dl.next_to(dl_elem, np.array([0,-7,0]))
        E1.move_to((-4.5, -2.5, 0))
        
        self.play(FadeIn(dim1_text), FadeIn(dom1), FadeIn(lamb), FadeIn(dl), FadeIn(dl_elem))
        
        
        # 2D
        
        dom2 = Rectangle(height = 3, width = 2, fill_opacity = opacity, color=color)
        dom2.next_to(dim2_text, np.array([0,-3,0]))
        dA_elem = Rectangle(height = 0.2, width = 0.2, fill_opacity = 1, color='white')
        dA_elem.move_to(dom2)
        E2 = TexMobject("E(\\boldsymbol{r}) = \\dfrac{1}{4\\pi \\epsilon_0} \\iint_{S'}  " + \
                             "\\dfrac{\\sigma dA}" +\
                                    "{\\boldsymbol{r}_s^2} \\hat{\\boldsymbol{r}}_s}")    
        
            
        sigma = TexMobject("dA")
        sigma.next_to(dA_elem)
        dA = TexMobject("dq = \\sigma dA")
        dA.next_to(dA_elem, np.array([0,-7,0]))
        E2.move_to((0, -2.5, 0))
        
        self.play(FadeIn(dim2_text), FadeIn(dom2), FadeIn(dA), FadeIn(dA_elem), FadeIn(sigma))
        
        # 3D
        
        dom3 = Ellipse(height = 3, width=2, fill_opacity = opacity, color=color)
        dom3.next_to(dim3_text, np.array([0,-3,0]))
        dV_elem = Circle(radius = 0.1 , fill_opacity = 1, color='white')
        dV_elem.move_to(dom3)
        E3 = TexMobject("E(\\boldsymbol{r}) = \\dfrac{1}{4\\pi \\epsilon_0} \\iiint_{V'}  " + \
                             "\\dfrac{\\rho dV}" +\
                                    "{\\boldsymbol{r}_s^2} \\hat{\\boldsymbol{r}}_s}")    
       
            
        rho = TexMobject("dV")
        rho.next_to(dV_elem)
        dV = TexMobject("dq = \\rho dV")
        dV.next_to(dV_elem, np.array([0,-7,0]))
        E3.next_to(E2)
        
        self.play(FadeIn(dim3_text),FadeIn(dom3), FadeIn(dV_elem), FadeIn(rho), FadeIn(dV))
        
        Epos = TexMobject(".")        
        Epos.move_to((3.5,3,0))
        
        E1.next_to(Epos, np.array([0, -2, 0]))
        E2.next_to(Epos, np.array([0,-10,0]))
        E3.next_to(Epos, np.array([0,-18,0]))
        
        self.play(FadeIn(E1), FadeIn(E2), FadeIn(E3))
        

class E_example(GraphScene):
    
    # 5. Example field calculation
    
    CONFIG = {
        "x_min": 0,
        "x_max": 1,
        "num_graph_anchor_points": 100,
        "y_min": -1,
        "y_max": 1,
        "graph_origin": np.array([-5,1,0]),
        "function_color": RED,
        "axes_color": GRAY,
        "x_axis_width": 5,
        "y_axis_width": 4.5
        #"x_labeled_nums": range(-10, 12, 2),
    }


    def construct(self):
        
        y_mid = np.array([-5,1,0])
        
        self.setup_axes()
        lineQ = Rectangle(height = 5, width = 0.2, fill_opacity = 0.5, color='cyan')
        lineQ.move_to(y_mid)
        
        P = Dot(color=GREEN)
        P_pos = (-0.5,1,0)
        P.move_to(P_pos)
        P_text = TextMobject("P")
        P_text.next_to(P, UP)
        
        dl_elem = Rectangle(height = 0.2, width = 0.2, fill_opacity = 1, color='white')
        dl_pos = np.array([-5,2.7,0])
        dl_elem.move_to(dl_pos)

        dyT = TextMobject("dy")
        dQT = TextMobject("dQ")
        
        dyT.next_to(dl_elem, LEFT)
        dQT.next_to(dl_elem, RIGHT)
        
        r_vec = Arrow(dl_pos, P_pos, stroke_width = 2)
        rV_text = TexMobject("\\boldsymbol{r}_s", height=0.4)
        rV_text.next_to(np.array([-2.5,1,0]), np.array([0,1,0]))
        
        rs_formula = TexMobject("\\boldsymbol{r}_s= x \\boldsymbol{i} - y\\boldsymbol{j}")
        rs_formula.move_to((3,2,0))
        
        self.play(FadeIn(lineQ), FadeIn(P), FadeIn(P_text))
        self.play(FadeIn(dl_elem),
                  FadeIn(dyT), FadeIn(dQT), FadeIn(r_vec), FadeIn(rV_text),
                  FadeIn(rs_formula))
        
        dfrac_size = 0.9
        E_formula = TexMobject("E(\\boldsymbol{r})= \\dfrac{1}{4\\pi \\epsilon_0}" + 
                               "\\int_{l'} \\dfrac{dQ}{r^2_s} \\boldsymbol{\\hat{r}}_s", height = dfrac_size)
        E_formula.next_to(rs_formula, DOWN)
        E2 = TexMobject("= \\dfrac{1}{4\\pi \\epsilon_0} \\int_{-a}^{a}"+
                        "\\dfrac{x \\boldsymbol{i} - y \\boldsymbol{j}}{(x^2+y^2)^{3/2}}"+
                        "\\lambda dy", height = dfrac_size)
        E2.move_to((-1,-0.5,0))
        
        E3 = TexMobject("E(\\boldsymbol{r})= \\dfrac{1}{4\\pi \\epsilon_0} \\left( \\dfrac{Q}{x\\sqrt{x^2+a^2}}, -0 \\right)" +
                        "= \\dfrac{1}{4\\pi \\epsilon_0} \\dfrac{Q}{x\\sqrt{x^2+a^2}} \\boldsymbol{i}", height=dfrac_size)
        E3.move_to((0,-2,0))
        ApplyMethod(E3.scale, 0.6)
        
        self.play(FadeIn(E_formula), FadeIn(E2))
        self.play(FadeIn(E3))
        #self.play(E3.scale, 0.6)

class ElecFlux(Scene):
    
    # 6. Electric flux definition.
    
    def construct(self):
        
        # Construct origin
        
        O_text = TextMobject("Sähkövuo")
        O_text.move_to(point_or_mobject = [0,3,0])
        self.play(FadeIn(O_text))   
        
        Q_center = np.array([-4, 0, 0])
        
        dom3 = Ellipse(height = 2.5, width=2, fill_opacity = 0.3, color='green')
        dom3.move_to(Q_center)
        Q_text = TexMobject("Q")
        Q_text.move_to(Q_center + np.array([0,1.5,0]))
        
        S_text = TexMobject("S")
        S_text.move_to(Q_center + np.array([-2.5,0,0]))
        surface = DashedVMobject(Circle(radius = 2))
        surface.move_to(Q_center)
        
        dS_pos = Q_center + np.array([2,0,0])
        dS_elem = Ellipse(height = 0.5, width = 0.2, fill_opacity=1, color='orange')
        dS_elem.move_to(dS_pos)
        
        dS_text = TexMobject("dS")
        dS_text.next_to(dS_elem)
        
        dS_vec = Arrow(dS_pos, dS_pos+np.array([1.5,0,0]))
        dS_text.next_to(dS_vec, UP)
        
        self.play(FadeIn(dom3), FadeIn(Q_text), ShowCreation(surface), 
                  FadeIn(dS_elem), FadeIn(dS_text), FadeIn(dS_vec),
                  FadeIn(S_text))
        
        # Formulas
        dfrac_size = 0.7
        text_start = np.array([3,3,0])
        flux_formula = TexMobject("\\Phi = \oint_S \\overline{E} \\cdot d \\overline{S}",
                                  height=dfrac_size)
        flux_formula.move_to(text_start)

        flux_formula2 = TexMobject(" = \oint_{S'} \\dfrac{1}{4\\pi\\epsilon_0} \\left[ \int_{V'} \\rho " +
                                   "\\dfrac{\\boldsymbol{r}-\\boldsymbol{r}'}{|\\boldsymbol{r}-\\boldsymbol{r}'|^3} dV \\right] \\cdot d \\overline{S}",
                                  height=dfrac_size)
        flux_formula2.move_to(text_start + np.array([0,-1.5,0]))
        
        flux_formula3 = TexMobject(" = \\dfrac{1}{4\\pi\\epsilon_0} \int_{V'} \\rho \\left[ \oint_{S'} " +
                                   "\\dfrac{\\boldsymbol{r}-\\boldsymbol{r}'}{|\\boldsymbol{r}-\\boldsymbol{r}'|^3}  \\cdot d \\overline{S} \\right] dV ",
                                  height=dfrac_size)
        flux_formula3.move_to(text_start + np.array([0,-2.5,0]))
        
        flux_formula4 = TexMobject(" = \\dfrac{1}{4\\pi\\epsilon_0} \int_{V'} \\rho 4 \\pi dV ",
                                  height=dfrac_size)
        flux_formula4.move_to(text_start + np.array([0,-3.5,0]))
        
        
        flux_formula5 = TexMobject(" = \\dfrac{1}{\\epsilon_0} \int_{V'} \\rho dV ",
                                  height=dfrac_size)
        flux_formula5.move_to(text_start + np.array([-0.25,-4.5,0]))
        
        flux_formula6 = TexMobject(" = \\dfrac{Q_{\\mathrm{enc}}}{\\epsilon_0} ",
                                  height=dfrac_size)
        final_pos = text_start + np.array([-0.5,-5.5,0])
        flux_formula6.move_to(final_pos)
        
        self.play(FadeIn(flux_formula))
        self.play(FadeIn(flux_formula2))
        self.play(FadeIn(flux_formula3))
        self.play(FadeIn(flux_formula4))
        self.play(FadeIn(flux_formula5))
        self.play(FadeIn(flux_formula6))
        #self.wait()
        
        self.play(FadeOut(flux_formula2), FadeOut(flux_formula3), FadeOut(flux_formula4),
                  FadeOut(flux_formula5))
        
        compare_pos = np.array([2,0,0])
        self.play(flux_formula.shift, compare_pos - text_start,
                  flux_formula6.shift, compare_pos - final_pos+np.array([1.6,0.04,0]))
        
        framebox1 = Rectangle(height = 2, width = 3.5,opacity=0, color = YELLOW)
        framebox1.move_to(compare_pos + np.array([0.55,0,0]))
        self.play(ShowCreation(framebox1))
        self.wait()
        
        gauss_text = TextMobject("Gaussin laki")
        gauss_text.next_to(framebox1, np.array([0,-1.5,0]))
        self.play(FadeIn(gauss_text))
         
class SolidAngle(Scene):
    
    def construct(self):
        
        # Construct origin
        
        O_text = TextMobject("Avaruuskulma", height = 0.5)
        O_text.move_to(point_or_mobject = [0,3,0])
        
        text_start = np.array([-3,1,0])
        
        sol_int = TexMobject("\oint_{S}\\dfrac{\\boldsymbol{r}-\\boldsymbol{r}'}"+
                             "{|\\boldsymbol{r}-\\boldsymbol{r}'|^3}  \\cdot d \\overline{S}"+
                             "=\oint_{S}\\dfrac{\\boldsymbol{\\hat{r}}_s \\cdot \\boldsymbol{\\hat{n}}}{r_s^2} dS",
                             height = 1)
        
        sol_int.move_to(text_start + np.array([0,0,0]))
        
        a = np.array([0,0,0])
        eq1A = TexMobject("4\\pi \\, \\, , \\, \\, \\overline{r}' \\in"
                          +"\\,V")
        eq2A = TexMobject("0 \\, \\, \\, , \\, \\, \\overline{r}' \\notin"
                          +"\\,V")
        
        #eq1A = TexMobject("4x + 3y")
        #eq2A = TexMobject("5x -2y")
       # eq1A.move_to((text_start + np.array([2,0,0])))
        eq1A.move_to(text_start + np.array([5.5,0.25,0]))
        #eq1A.next_to(sol_int)
        eq2A.next_to(eq1A,np.array([0,-1,0]))
        eqE = TexMobject("=")
        #eq1C.next_to(eq1B,RIGHT)
        #eq2A.shift(DOWN)
        #eq2B.shift(DOWN)
        #eq2C.shift(DOWN)
        #eq2A.align_to(eq1A,LEFT)
        #eq2B.align_to(eq1B,LEFT)
        #eq2C.align_to(eq1C,LEFT)

        eq_group=VGroup(eq1A,eq2A)
        braces=Brace(eq_group,LEFT)
        eqE.next_to(braces, LEFT)
        #eq_text = braces.get_text("A pair of equations")

        self.play(FadeIn(O_text))
        self.play(FadeIn(sol_int))
        self.play(FadeIn(eq1A), FadeIn(eq2A), FadeIn(eqE), GrowFromCenter(braces))
        #self.play(GrowFromCenter(braces))
        
        #
        
        r1_pos = np.array([-5,-3,0])
        r1d_pos = np.array([-2,-2,0,])
        r2_pos = np.array([0,-3,0])
        r2d_pos = np.array([3,-2,0,])        
        
        vec1 = Arrow(r1_pos, r1d_pos, stroke_width = 2)
        vec2 = Arrow(r2_pos, r2d_pos, stroke_width = 2)
        rs1 = TexMobject("\\boldsymbol{r}'")
        rs1.move_to(r1d_pos - np.array([1.5,0,0]))
        rs2 = TexMobject("\\boldsymbol{r}'")
        rs2.move_to(r2d_pos - np.array([1.5,0,0]))
        
        S1 = DashedVMobject(Circle(radius = 0.8))
        S1.move_to(r1d_pos)
        S2 = DashedVMobject(Circle(radius = 0.8))
        S2.move_to(r2d_pos + np.array([1,0,0]))
        
        zero = TexMobject("0")
        fPi = TexMobject("4\\pi")
        zero.move_to(r2d_pos + np.array([1, 1.5,0]))
        fPi.move_to(r1d_pos + np.array([0,1.5,0]))
        
        self.play(FadeIn(vec1), FadeIn(vec2), FadeIn(rs1), FadeIn(rs2),
                  ShowCreation(S1), ShowCreation(S2), FadeIn(zero), FadeIn(fPi))
        
        
        
class GaussExample(ThreeDScene):
    
    # 8. Gauss law example
    
    def construct(self):
        resolution_fa = 4
        self.set_camera_orientation(phi=75 * DEGREES, theta=-30 * DEGREES)

        def param_plane(u, v):
            x = u
            y = v
            z = 0
            return np.array([x, y, z])

        max_values = [-2,2,-2,2]
        min_u = -2
        min_v = -2
        max_u = 2
        max_v = 2

        plane = ParametricSurface(
            param_plane,
            resolution=(resolution_fa, resolution_fa),
            v_min= min_v,
            v_max= max_v,
            u_min= min_u,
            u_max= max_u,
        )
        plane.set_fill_by_checkerboard(YELLOW, opacity=0.5)
        
        axes = ThreeDAxes()
        scale = 1.5
        plane.scale_about_point(scale, ORIGIN)

        
        #myVec=CurvedArrow(start_point=np.array([2,2,0]),end_point=np.array([2,2,2]),color=GOLD_D)
        
        
        self.set_camera_orientation(phi=55 * DEGREES, theta=35 * DEGREES)
        self.add(axes)
        self.play(FadeIn(plane))
        
        for i in range(0, resolution_fa) :
            for j in range(0,resolution_fa):
                plus = TextMobject("+", height = 0.4)
                # min_u + 2/num_rect + i*(max_u-min_u)/num_rect
                plus.move_to(np.array([scale*(min_u + 2/resolution_fa + i*(max_u-min_u)/resolution_fa)
                                       ,scale*(min_v + 2/resolution_fa + j*(max_v-min_v)/resolution_fa)
                                       , 0]))
                self.add(plus)
        
        self.begin_ambient_camera_rotation(rate=0.3)
        #self.play(Rotate(arrow, -PI/2, axis=Y_AXIS))
        self.wait(6)
        self.stop_ambient_camera_rotation()
        #self.move_camera(phi=90 * DEGREES, theta=60 * DEGREES)
        arrow1 = Arrow(np.array([1,2,3]),  np.array([3,2,3]), color = WHITE)
        arrow2 = Arrow(np.array([1,2,-3]),  np.array([3,2,-3]), color = WHITE)
        arrow1.rotate(-PI/2, axis = Y_AXIS) 
        arrow2.rotate(PI/2, axis = Y_AXIS) 

        
        # Gaussian cylinder
        
        def param_cylinder(u, v):
            x = np.cos(u)
            y = np.sin(u)
            z = v
            return np.array([x, y, z])

        cylinder = ParametricSurface(
            param_cylinder,
            resolution=(10, 20),
            v_min= -0.5,
            v_max= 0.5,
            u_min= 0,
            u_max= 2*np.pi,
        )
        
        cylinder_pos = np.array([2.6,3.2,-1])
        top_circle = Circle(radius = 1, color = GREEN, fill_opacity=0.5)
        top_circle.move_to(cylinder_pos + np.array([0,0, 0.5]))
        bot_circle = Circle(radius = 1, color = GREEN, fill_opacity=0.5)
        bot_circle.move_to(cylinder_pos + np.array([0,0, 0.5]))
        
        cylinder.set_fill_by_checkerboard(BLUE, opacity=0.8)
        cylinder.move_to(cylinder_pos)
        
        #self.move_camera(phi=75 * DEGREES, theta=60 * DEGREES)
        
        arrow1.move_to(cylinder_pos + np.array([0,0,1+0.5]))
        arrow2.move_to(cylinder_pos + np.array([0,0,-1-0.5]))
        self.wait(2)
        self.play(FadeIn(arrow2),
                  FadeIn(cylinder), FadeIn(bot_circle), FadeIn(top_circle),
                  FadeIn(arrow1))
        

        text_pos = np.array([3,3,0])
        E_text1 = TexMobject("\\int_S \\boldsymbol{E} \\cdot \\hat{\\boldsymbol{z}} dA = \\dfrac{Q_{\\mathrm{enc}}}{\\epsilon_0} ")
        E_text2 = TexMobject(" = \\dfrac{\\sigma A}{\\epsilon_0}") # \\dfrac{\\sigma A}{\\epsilon_0}
        E_text3 = TexMobject("2 A |\\boldsymbol{E}| ")
        E_text4 = TexMobject("\\boldsymbol{E}=\\dfrac{\\sigma}{2\\epsilon_0} \\hat{\\boldsymbol{z}}")
        
        E_text1.move_to(text_pos)
        E_text2.move_to(np.array([4.5,1.7,0]))
        E_text3.move_to(np.array([4.5,0.5,0]))
        E_text4.move_to(np.array([4.5,-1,0]))
        self.add_fixed_in_frame_mobjects(E_text1)
        self.play(FadeIn(E_text1) )
        self.add_fixed_in_frame_mobjects(E_text2)
        self.play(FadeIn(E_text2))
        self.add_fixed_in_frame_mobjects(E_text3)
        self.play(FadeIn(E_text3))
        self.add_fixed_in_frame_mobjects(E_text4)
        self.play(FadeIn(E_text4))
        #self.play(FadeIn(E_text2))

        
class CylindricalCoordinates(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=45 * DEGREES)
        self.camera.frame_center.shift(3.5 * LEFT)
        axes = ThreeDAxes()
        cylinder = ParametricSurface(
            lambda u, v: np.array([
                np.cos(TAU * v),
                np.sin(TAU * v),
                2 * (1 - u)
            ]),
            resolution=(6, 32)).fade(0.5)  # Resolution of the surfaces
        self.add(cylinder.scale(2))
        dot = Dot()
        def func(t):
            return np.array((np.cos(2* t), np.sin(2* t), 0.5+np.cos(t)))
        func = ParametricFunction(func, t_max=TAU, fill_opacity=0)
        func.scale(2)
        new_func = CurvesAsSubmobjects(func.copy())
        new_func.set_color_by_gradient(BLUE, WHITE, RED)
        dot.add_updater(lambda m: m.move_to(func.get_end()))
        func.fade(1)
        self.add(dot)
        text = TexMobject(r" S_1(",r"t",r")= R \left( \begin{array}{c}  \cos(2 \omega t) \\ \sin(2 \omega t)\\ 0 \end{array} \right)+ \left( \begin{array}{c} 0 \\ 0\\ \cos(\omega t) \end{array} \right) ")
        text.to_corner(LEFT).shift(3.5*LEFT).scale(0.8)
        text2 = text.copy()
        text[1].set_color(BLUE)
        self.add_fixed_in_frame_mobjects(text)
        text2[1].set_color(RED)
        text3= TexMobject("t = [","0",",","  \\frac{2\pi}{\omega}","]")
        text3[1].set_color(BLUE)
        text3[3].set_color(RED)
        text3.next_to(text2,DOWN)
        self.add_fixed_in_frame_mobjects(text3)
        self.wait()
        self.play(
            ShowCreation(new_func),
            ShowCreation(func),
            Transform(text, text2),
            run_time=3.5
        )
        self.wait()
        
class EFieldInThreeD(ThreeDScene):
    CONFIG = {
        "plane_kwargs" : {
            "color" : RED_B
            },
        "point_charge_loc" : 0.5*RIGHT-1.5*UP,
        }
    def calc_field2D(self,point):
        x,y = point[:2]
        Rx,Ry = self.point_charge_loc[:2]
        r = math.sqrt((x-Rx)**2 + (y-Ry)**2)
        efield = (point - self.point_charge_loc)/r**3
        return Vector(efield).shift(point)
 
    def calc_field3D(self,point):
        x,y,z = point
        Rx,Ry,Rz = self.point_charge_loc
        r = math.sqrt((x-Rx)**2 + (y-Ry)**2+(z-Rz)**2)
        efield = (point - self.point_charge_loc)/r**3
        return Vector(efield).shift(point)
        
    
    def construct(self):
        
        self.set_camera_orientation(phi=0 * DEGREES, theta=-30 * DEGREES)
        
        O = TexMobject("lol")
        self.play(FadeIn(O))
        
        testi = self.calc_field3D(np.array([2,2,2]))
        testi2 = Vector(np.array([2,2,2]))
        
        self.play(FadeIn(testi))
        #self.play(FadeIn(testi2))
        
        vector=Vector(direction=0.5*(UP+RIGHT+OUT))
        self.add(vector)
        self.wait(2)
 


class VectorFieldScene1(Scene):
    def construct(self):
        func = lambda p: np.array([
            p[0]/2,  # x
            p[1]/2,  # y
            0        # z
        ])
        # Normalized
        vector_field_norm = VectorField(func)
        # Not normalized
        vector_field_not_norm = VectorField(func, length_func=linear)
        self.play(*[GrowArrow(vec) for vec in vector_field_norm])
        self.wait(2)
        self.play(ReplacementTransform(vector_field_norm,vector_field_not_norm))
        self.wait(2)
        
class TestScene(ThreeDScene):
    def construct(self):
        self.wait()
        prompt=TextMobject("Text here").to_corner(UP)
        cube = VGroup(
            Line(np.array([-1,-1,-1]),np.array([1,-1,-1])),
            Line(np.array([1,-1,-1]),np.array([1,-1,1])),
            Line(np.array([1,-1,1]),np.array([-1,-1,1])),
            Line(np.array([-1,-1,1]),np.array([-1,-1,-1])),
            Line(np.array([-1,1,-1]),np.array([1,1,-1])),
            Line(np.array([1,1,-1]),np.array([1,1,1])),
            Line(np.array([1,1,1]),np.array([-1,1,1])),
            Line(np.array([-1,1,1]),np.array([-1,1,-1])),
            Line(np.array([-1,-1,-1]),np.array([-1,1,-1])),
            Line(np.array([1,1,-1]),np.array([1,-1,-1])),
            Line(np.array([1,-1, 1]),np.array([1,1,1])),
            Line(np.array([-1,-1,1]),np.array([-1,1,1])),
            )
        self.set_camera_orientation(0.8*np.pi/2, -0.45*np.pi)
        self.add_fixed_in_frame_mobjects(prompt) # <---- Add this line
        self.play(
            LaggedStartMap(ShowCreation,cube.scale(2)),
            )
        self.wait()
        self.begin_ambient_camera_rotation(rate=0.3)
        self.wait(4)
        self.play(Write(prompt))
        