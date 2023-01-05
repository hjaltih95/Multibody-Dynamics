function Acc = acchw7(alpha,alphad,betad,beta,gammad,gamma)
%ACCHW7
%    ACC = ACCHW7(ALPHA,ALPHAD,BETAD,BETA,GAMMAD,GAMMA)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Jul-2019 21:30:23

t2 = sin(alpha);
t3 = cos(gamma);
t4 = t3.^2;
t5 = cos(beta);
t6 = sin(beta);
t7 = t6.^2;
t8 = t2.^2;
t9 = t7.^2;
t10 = t5.^2;
t11 = sin(gamma);
t12 = t11.^2;
t13 = beta-gamma;
t14 = cos(t13);
t15 = beta+gamma;
t16 = cos(t15);
t17 = betad.^2;
t18 = alphad.^2;
t19 = cos(alpha);
t20 = gammad.^2;
t21 = t4.*1.215243e6;
t22 = t2.*t4.*3.240648e6;
t23 = t9.*t10.*6.784e3;
t24 = t4.*t8.*3.240648e6;
t25 = t2.*t4.*t8.*1.440288e6;
t26 = t8.^2;
t27 = t4.*t26.*2.40048e5;
t28 = t4.*t9.*3.1008e4;
t29 = t4.*t7.*t9.*2.352e3;
t30 = t2.*t9.*t10.*2.0352e4;
t31 = t2.*t4.*t9.*9.3024e4;
t32 = t2.*t4.*t7.*t9.*7.056e3;
t33 = t8.*t9.*t10.*1.5264e4;
t34 = t4.*t8.*t9.*6.9768e4;
t35 = t4.*t7.*t8.*t9.*5.292e3;
t36 = t9.*t10.*t12.*7.68e2;
t37 = t8.*t9.*t10.*t12.*1.728e3;
t38 = t3.*t5.*t9.*t16.*1.344e3;
t39 = t2.*t9.*t10.*t12.*2.304e3;
t40 = t2.*t3.*t5.*t9.*t16.*4.032e3;
t41 = t3.*t5.*t8.*t9.*t16.*3.024e3;
t45 = t3.*t5.*t9.*t14.*4.032e3;
t46 = t3.*t5.*t8.*t9.*t14.*9.072e3;
t47 = t2.*t3.*t5.*t9.*t14.*1.2096e4;
t42 = t21+t22+t23+t24+t25+t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41-t45-t46-t47;
t43 = 1.0./t42;
t44 = 1.0./t6;
t48 = t9.^2;
t49 = t16.^2;
t50 = t14.^2;
t51 = t2.*2.0;
t52 = t51+3.0;
t53 = t5.*t7.*t17.*1.1448e4;
t54 = t5.*t7.*t8.*t17.*1.5264e4;
t55 = t2.*t5.*t7.*t8.*t17.*3.392e3;
t56 = t5.*t7.*t12.*t17.*1.296e3;
t57 = t5.*t6.*t11.*t20.*1.7172e4;
t58 = t3.*t5.*t6.*t11.*1.62e6;
t59 = t3.*t7.*t16.*t17.*1.134e3;
t60 = t2.*t5.*t7.*t17.*2.2896e4;
t61 = t5.*t6.*t11.*t12.*t20.*1.944e3;
t62 = t3.*t10.*t14.*t17.*1.3608e4;
t63 = t2.*t6.*t7.*t10.*t18.*t19.*5.088e3;
t64 = t2.*t4.*t6.*t7.*t18.*t19.*2.3256e4;
t65 = t2.*t4.*t6.*t9.*t18.*t19.*1.764e3;
t66 = t2.*t5.*t7.*t12.*t17.*2.592e3;
t67 = t3.*t6.*t11.*t16.*t20.*1.701e3;
t68 = t5.*t6.*t8.*t11.*t12.*t20.*2.592e3;
t69 = t2.*t5.*t6.*t8.*t11.*t12.*t20.*5.76e2;
t70 = t3.*t8.*t10.*t14.*t17.*1.8144e4;
t71 = t2.*t3.*t8.*t10.*t14.*t17.*4.032e3;
t72 = t2.*t5.*t6.*t11.*t20.*3.4344e4;
t73 = t2.*t3.*t5.*t6.*t11.*3.24e6;
t74 = t2.*t3.*t7.*t16.*t17.*2.268e3;
t75 = t5.*t7.*t8.*t12.*t17.*1.728e3;
t76 = t2.*t5.*t7.*t8.*t12.*t17.*3.84e2;
t77 = t3.*t6.*t10.*t11.*t17.*4.536e3;
t78 = t4.*t5.*t6.*t11.*t20.*5.346e3;
t79 = t5.*t6.*t8.*t11.*t20.*2.2896e4;
t80 = t2.*t5.*t6.*t11.*t12.*t20.*3.888e3;
t81 = t2.*t5.*t6.*t8.*t11.*t20.*5.088e3;
t82 = t2.*t3.*t10.*t14.*t17.*2.7216e4;
t83 = t3.*t5.*t6.*t8.*t11.*2.16e6;
t84 = t2.*t3.*t5.*t6.*t8.*t11.*4.8e5;
t85 = t3.*t7.*t8.*t16.*t17.*1.512e3;
t86 = t2.*t3.*t7.*t8.*t16.*t17.*3.36e2;
t87 = t2.*t3.*t6.*t11.*t16.*t20.*3.402e3;
t88 = t2.*t6.*t7.*t10.*t12.*t18.*t19.*5.76e2;
t89 = t3.*t6.*t8.*t11.*t16.*t20.*2.268e3;
t90 = t2.*t3.*t6.*t8.*t11.*t16.*t20.*5.04e2;
t91 = t2.*t3.*t6.*t10.*t11.*t17.*9.072e3;
t92 = t3.*t5.*t6.*t7.*t14.*t18.*t19.*5.04e2;
t93 = t2.*t4.*t5.*t6.*t11.*t20.*1.0692e4;
t94 = t3.*t6.*t8.*t10.*t11.*t17.*6.048e3;
t95 = t2.*t3.*t6.*t8.*t10.*t11.*t17.*1.344e3;
t96 = t4.*t5.*t6.*t8.*t11.*t20.*7.128e3;
t97 = t2.*t4.*t5.*t6.*t8.*t11.*t20.*1.584e3;
t98 = t2.*t3.*t5.*t6.*t7.*t16.*t18.*t19.*1.008e3;
t99 = t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98-t4.*t7.*2.835e6-t2.*t4.*t7.*5.67e6-t4.*t7.*t8.*3.78e6-t4.*t5.*t17.*1.04652e5-t5.*t10.*t17.*2.2896e4-t2.*t4.*t7.*t8.*8.4e5-t2.*t4.*t5.*t17.*2.09304e5-t4.*t5.*t7.*t17.*1.4175e4-t2.*t5.*t10.*t17.*4.5792e4-t3.*t4.*t7.*t20.*6.804e3-t4.*t5.*t8.*t17.*1.39536e5-t5.*t8.*t10.*t17.*3.0528e4-t3.*t7.*t14.*t17.*3.402e3-t5.*t10.*t12.*t17.*2.592e3-t3.*t10.*t16.*t17.*4.536e3-t2.*t4.*t5.*t7.*t17.*2.835e4-t2.*t3.*t4.*t7.*t20.*1.3608e4-t2.*t4.*t5.*t8.*t17.*3.1008e4-t4.*t5.*t7.*t8.*t17.*1.89e4-t2.*t5.*t8.*t10.*t17.*6.784e3-t3.*t4.*t7.*t8.*t20.*9.072e3-t2.*t3.*t7.*t14.*t17.*6.804e3-t2.*t5.*t10.*t12.*t17.*5.184e3-t2.*t3.*t10.*t16.*t17.*9.072e3-t3.*t7.*t8.*t14.*t17.*4.536e3-t5.*t8.*t10.*t12.*t17.*3.456e3-t3.*t6.*t11.*t14.*t20.*5.103e3-t3.*t8.*t10.*t16.*t17.*6.048e3-t4.*t6.*t7.*t18.*t19.*3.876e3-t4.*t6.*t9.*t18.*t19.*2.94e2-t6.*t7.*t10.*t18.*t19.*8.48e2-t2.*t4.*t5.*t7.*t8.*t17.*4.2e3-t2.*t3.*t4.*t7.*t8.*t20.*2.016e3-t2.*t3.*t7.*t8.*t14.*t17.*1.008e3-t2.*t5.*t8.*t10.*t12.*t17.*7.68e2-t2.*t3.*t6.*t11.*t14.*t20.*1.0206e4-t2.*t3.*t8.*t10.*t16.*t17.*1.344e3-t3.*t6.*t8.*t11.*t14.*t20.*6.804e3-t6.*t7.*t10.*t12.*t18.*t19.*9.6e1-t2.*t3.*t6.*t8.*t11.*t14.*t20.*1.512e3-t3.*t5.*t6.*t7.*t16.*t18.*t19.*1.68e2-t2.*t3.*t5.*t6.*t7.*t14.*t18.*t19.*3.024e3;
Acc = [(t6.*t43.*t99.*(t2.*3.0+2.0).*2.0)./t52;t43.*t44.*(t4.*t7.*t9.*-8.4e5+t4.*t5.*t17.*1.215243e6-t2.*t4.*t7.*t9.*2.52e6+t2.*t4.*t5.*t17.*3.240648e6-t4.*t7.*t8.*t9.*1.89e6+t4.*t5.*t8.*t17.*3.240648e6+t5.*t7.*t9.*t17.*3.392e3+t4.*t5.*t17.*t26.*2.40048e5+t3.*t5.*t6.*t9.*t11.*4.8e5+t2.*t4.*t5.*t8.*t17.*1.440288e6+t2.*t5.*t7.*t9.*t17.*1.0176e4-t4.*t5.*t7.*t9.*t17.*1.848e3-t3.*t4.*t7.*t9.*t20.*2.016e3+t5.*t7.*t8.*t9.*t17.*7.632e3-t3.*t7.*t9.*t14.*t17.*1.008e3+t5.*t7.*t9.*t12.*t17.*3.84e2+t5.*t6.*t9.*t11.*t20.*5.088e3+t3.*t7.*t9.*t16.*t17.*3.36e2+t4.*t6.*t7.*t18.*t19.*4.5009e4+t2.*t3.*t5.*t6.*t9.*t11.*1.44e6+t3.*t5.*t6.*t8.*t9.*t11.*1.08e6-t2.*t4.*t5.*t7.*t9.*t17.*5.544e3-t2.*t3.*t4.*t7.*t9.*t20.*6.048e3-t4.*t5.*t7.*t8.*t9.*t17.*4.158e3-t3.*t4.*t7.*t8.*t9.*t20.*4.536e3-t2.*t3.*t7.*t9.*t14.*t17.*3.024e3+t2.*t5.*t7.*t9.*t12.*t17.*1.152e3+t2.*t5.*t6.*t9.*t11.*t20.*1.5264e4+t2.*t3.*t7.*t9.*t16.*t17.*1.008e3+t4.*t5.*t6.*t9.*t11.*t20.*1.584e3-t2.*t4.*t6.*t7.*t18.*t19.*2.40048e5+t3.*t6.*t9.*t10.*t11.*t17.*1.344e3-t3.*t7.*t8.*t9.*t14.*t17.*2.268e3+t5.*t7.*t8.*t9.*t12.*t17.*8.64e2+t5.*t6.*t8.*t9.*t11.*t20.*1.1448e4+t3.*t7.*t8.*t9.*t16.*t17.*7.56e2-t4.*t6.*t7.*t8.*t18.*t19.*1.80036e5-t3.*t6.*t9.*t11.*t14.*t20.*1.512e3+t5.*t6.*t9.*t11.*t12.*t20.*5.76e2+t3.*t6.*t9.*t11.*t16.*t20.*5.04e2+t2.*t4.*t5.*t6.*t9.*t11.*t20.*4.752e3+t2.*t3.*t6.*t9.*t10.*t11.*t17.*4.032e3+t4.*t5.*t6.*t8.*t9.*t11.*t20.*3.564e3+t3.*t6.*t8.*t9.*t10.*t11.*t17.*3.024e3-t2.*t3.*t6.*t9.*t11.*t14.*t20.*4.536e3+t2.*t5.*t6.*t9.*t11.*t12.*t20.*1.728e3+t2.*t3.*t6.*t9.*t11.*t16.*t20.*1.512e3-t3.*t6.*t8.*t9.*t11.*t14.*t20.*3.402e3+t5.*t6.*t8.*t9.*t11.*t12.*t20.*1.296e3+t3.*t6.*t8.*t9.*t11.*t16.*t20.*1.134e3).*2.0;t43.*t44.*(t3.*t7.*t17.*8.10162e5-t3.*t10.*t17.*1.620324e6+t3.*t17.*t48.*1.568e3+t3.*t5.*t7.*t9.*1.12e6+t2.*t3.*t7.*t17.*2.160432e6-t2.*t3.*t10.*t17.*4.320864e6+t3.*t7.*t8.*t17.*2.160432e6+t3.*t7.*t9.*t17.*2.0672e4-t6.*t9.*t10.*t11.*6.4e5-t3.*t8.*t10.*t17.*4.320864e6+t3.*t6.*t11.*t20.*1.215243e6+t3.*t7.*t17.*t26.*1.60032e5-t3.*t10.*t17.*t26.*3.20064e5+t2.*t3.*t17.*t48.*4.704e3+t3.*t8.*t17.*t48.*3.528e3+t2.*t3.*t5.*t7.*t9.*3.36e6+t3.*t5.*t7.*t8.*t9.*2.52e6+t2.*t3.*t7.*t8.*t17.*9.60192e5+t2.*t3.*t7.*t9.*t17.*6.2016e4-t2.*t6.*t9.*t10.*t11.*1.92e6-t2.*t3.*t8.*t10.*t17.*1.920384e6+t2.*t3.*t6.*t11.*t20.*3.240648e6+t3.*t7.*t8.*t9.*t17.*4.6512e4-t6.*t8.*t9.*t10.*t11.*1.44e6+t4.*t5.*t7.*t9.*t20.*2.688e3+t3.*t7.*t9.*t10.*t17.*2.464e3+t3.*t6.*t8.*t11.*t20.*3.240648e6+t3.*t6.*t9.*t11.*t20.*3.1008e4-t5.*t7.*t9.*t14.*t17.*1.344e3+t5.*t7.*t9.*t16.*t17.*4.48e2+t3.*t6.*t11.*t20.*t26.*2.40048e5+t2.*t4.*t5.*t7.*t9.*t20.*8.064e3+t2.*t3.*t7.*t9.*t10.*t17.*7.392e3+t2.*t3.*t6.*t8.*t11.*t20.*1.440288e6+t2.*t3.*t6.*t9.*t11.*t20.*9.3024e4+t4.*t5.*t7.*t8.*t9.*t20.*6.048e3-t2.*t5.*t7.*t9.*t14.*t17.*4.032e3+t3.*t7.*t8.*t9.*t10.*t17.*5.544e3+t2.*t5.*t7.*t9.*t16.*t17.*1.344e3+t3.*t6.*t7.*t9.*t11.*t20.*2.352e3+t3.*t6.*t8.*t9.*t11.*t20.*6.9768e4-t3.*t5.*t6.*t7.*t18.*t19.*6.0012e4-t5.*t6.*t9.*t10.*t11.*t17.*1.792e3-t3.*t6.*t9.*t10.*t11.*t20.*2.112e3-t5.*t7.*t8.*t9.*t14.*t17.*3.024e3+t5.*t7.*t8.*t9.*t16.*t17.*1.008e3-t5.*t6.*t9.*t11.*t14.*t20.*2.016e3+t5.*t6.*t9.*t11.*t16.*t20.*6.72e2+t2.*t3.*t6.*t7.*t9.*t11.*t20.*7.056e3+t2.*t3.*t5.*t6.*t7.*t18.*t19.*3.20064e5-t2.*t5.*t6.*t9.*t10.*t11.*t17.*5.376e3-t2.*t3.*t6.*t9.*t10.*t11.*t20.*6.336e3+t3.*t6.*t7.*t8.*t9.*t11.*t20.*5.292e3+t3.*t5.*t6.*t7.*t8.*t18.*t19.*2.40048e5-t5.*t6.*t8.*t9.*t10.*t11.*t17.*4.032e3-t2.*t5.*t6.*t9.*t11.*t14.*t20.*6.048e3-t3.*t6.*t8.*t9.*t10.*t11.*t20.*4.752e3+t2.*t5.*t6.*t9.*t11.*t16.*t20.*2.016e3-t5.*t6.*t8.*t9.*t11.*t14.*t20.*4.536e3+t5.*t6.*t8.*t9.*t11.*t16.*t20.*1.512e3);t43.*t44.*(t7.*t17.*4.2938586e7-t10.*t17.*8.5877172e7+t17.*t48.*8.3104e4+t3.*t6.*t11.*6.076215e9+t5.*t7.*t9.*5.936e7+t2.*t7.*t17.*1.14502896e8-t4.*t7.*t17.*8.506701e6-t2.*t10.*t17.*2.29005792e8+t7.*t8.*t17.*1.14502896e8+t7.*t9.*t17.*1.095616e6-t8.*t10.*t17.*2.29005792e8+t7.*t12.*t17.*4.860972e6+t6.*t11.*t20.*6.4407879e7-t10.*t12.*t17.*9.721944e6+t7.*t17.*t26.*8.481696e6-t10.*t17.*t26.*1.6963392e7+t2.*t17.*t48.*2.49312e5-t4.*t17.*t48.*1.6464e4+t8.*t17.*t48.*1.86984e5+t12.*t17.*t48.*9.408e3+t2.*t3.*t6.*t11.*1.620324e10+t2.*t5.*t7.*t9.*1.7808e8+t3.*t6.*t8.*t11.*1.620324e10+t3.*t6.*t9.*t11.*1.5504e8+t5.*t7.*t8.*t9.*1.3356e8-t2.*t4.*t7.*t17.*2.2684536e7-t3.*t7.*t9.*t14.*1.764e7+t5.*t7.*t9.*t12.*6.72e6+t2.*t7.*t8.*t17.*5.0890176e7+t2.*t7.*t9.*t17.*3.286848e6+t3.*t7.*t9.*t16.*5.88e6-t4.*t7.*t8.*t17.*2.2684536e7-t2.*t8.*t10.*t17.*1.01780352e8-t4.*t7.*t9.*t17.*2.17056e5+t2.*t7.*t12.*t17.*1.2962592e7+t2.*t6.*t11.*t20.*1.71754344e8+t3.*t5.*t14.*t17.*2.5520103e7-t2.*t10.*t12.*t17.*2.5925184e7-t3.*t5.*t16.*t17.*8.506701e6+t4.*t6.*t11.*t20.*7.291458e6+t7.*t8.*t9.*t17.*2.465136e6+t7.*t9.*t10.*t17.*8.3104e4+t7.*t8.*t12.*t17.*1.2962592e7+t6.*t8.*t11.*t20.*1.71754344e8+t7.*t9.*t12.*t17.*1.24032e5+t3.*t6.*t11.*t26.*1.20024e9+t6.*t9.*t11.*t20.*1.643424e6-t8.*t10.*t12.*t17.*2.5925184e7+t6.*t11.*t12.*t20.*7.291458e6-t4.*t7.*t17.*t26.*1.680336e6+t7.*t12.*t17.*t26.*9.60192e5+t6.*t11.*t20.*t26.*1.2722544e7-t10.*t12.*t17.*t26.*1.920384e6-t2.*t4.*t17.*t48.*4.9392e4-t4.*t8.*t17.*t48.*3.7044e4+t2.*t12.*t17.*t48.*2.8224e4-t7.*t9.*t17.*t49.*2.352e3-t7.*t9.*t17.*t50.*2.1168e4+t8.*t12.*t17.*t48.*2.1168e4+t2.*t3.*t6.*t8.*t11.*7.20144e9+t2.*t3.*t6.*t9.*t11.*4.6512e8-t2.*t3.*t7.*t9.*t14.*5.292e7+t2.*t5.*t7.*t9.*t12.*2.016e7+t3.*t6.*t7.*t9.*t11.*1.176e7+t2.*t3.*t7.*t9.*t16.*1.764e7+t3.*t6.*t8.*t9.*t11.*3.4884e8-t2.*t4.*t7.*t8.*t17.*1.0082016e7-t2.*t4.*t7.*t9.*t17.*6.51168e5+t2.*t3.*t5.*t14.*t17.*6.8053608e7-t3.*t7.*t8.*t9.*t14.*3.969e7+t5.*t7.*t8.*t9.*t12.*1.512e7+t3.*t5.*t6.*t11.*t17.*1.7013402e7-t2.*t3.*t5.*t16.*t17.*2.2684536e7+t2.*t4.*t6.*t11.*t20.*1.9443888e7+t3.*t7.*t8.*t9.*t16.*1.323e7+t3.*t5.*t7.*t9.*t20.*1.42464e5+t2.*t7.*t9.*t10.*t17.*2.49312e5-t4.*t7.*t8.*t9.*t17.*4.88376e5-t5.*t6.*t9.*t11.*t14.*1.008e7+t2.*t7.*t8.*t12.*t17.*5.761152e6+t2.*t6.*t8.*t11.*t20.*7.6335264e7+t2.*t7.*t9.*t12.*t17.*3.72096e5+t3.*t5.*t8.*t14.*t17.*6.8053608e7+t5.*t6.*t9.*t11.*t16.*3.36e6+t2.*t6.*t9.*t11.*t20.*4.930272e6-t2.*t8.*t10.*t12.*t17.*1.1522304e7-t3.*t5.*t8.*t16.*t17.*2.2684536e7+t4.*t6.*t8.*t11.*t20.*1.9443888e7+t4.*t6.*t9.*t11.*t20.*1.86048e5+t2.*t6.*t11.*t12.*t20.*1.9443888e7+t7.*t8.*t9.*t10.*t17.*1.86984e5+t6.*t7.*t9.*t11.*t20.*1.24656e5+t7.*t8.*t9.*t12.*t17.*2.79072e5-t4.*t7.*t9.*t14.*t20.*4.2336e4+t6.*t8.*t9.*t11.*t20.*3.697704e6-t5.*t6.*t7.*t18.*t19.*3.180636e6+t7.*t9.*t10.*t12.*t17.*9.408e3+t4.*t7.*t9.*t16.*t20.*1.4112e4-t6.*t9.*t10.*t11.*t20.*7.1232e4+t6.*t8.*t11.*t12.*t20.*1.9443888e7+t6.*t9.*t11.*t12.*t20.*1.86048e5+t7.*t9.*t14.*t16.*t17.*1.4112e4+t3.*t5.*t14.*t17.*t26.*5.041008e6-t3.*t5.*t16.*t17.*t26.*1.680336e6+t4.*t6.*t11.*t20.*t26.*1.440288e6+t6.*t11.*t12.*t20.*t26.*1.440288e6-t2.*t7.*t9.*t17.*t49.*7.056e3-t2.*t7.*t9.*t17.*t50.*6.3504e4-t7.*t8.*t9.*t17.*t49.*5.292e3-t7.*t8.*t9.*t17.*t50.*4.7628e4-t6.*t9.*t11.*t20.*t49.*3.528e3-t6.*t9.*t11.*t20.*t50.*3.1752e4+t2.*t3.*t6.*t7.*t9.*t11.*3.528e7+t2.*t3.*t5.*t6.*t11.*t17.*4.5369072e7+t3.*t6.*t7.*t8.*t9.*t11.*2.646e7+t2.*t3.*t5.*t7.*t9.*t20.*4.27392e5-t2.*t5.*t6.*t9.*t11.*t14.*3.024e7+t2.*t3.*t5.*t8.*t14.*t17.*3.0246048e7+t2.*t5.*t6.*t9.*t11.*t16.*1.008e7+t3.*t5.*t6.*t8.*t11.*t17.*4.5369072e7-t2.*t3.*t5.*t8.*t16.*t17.*1.0082016e7+t2.*t4.*t6.*t8.*t11.*t20.*8.641728e6+t3.*t5.*t6.*t9.*t11.*t17.*4.34112e5+t2.*t4.*t6.*t9.*t11.*t20.*5.58144e5+t3.*t5.*t7.*t8.*t9.*t20.*3.20544e5-t5.*t6.*t8.*t9.*t11.*t14.*2.268e7+t2.*t6.*t7.*t9.*t11.*t20.*3.73968e5-t3.*t5.*t7.*t9.*t14.*t17.*1.0584e4+t5.*t6.*t8.*t9.*t11.*t16.*7.56e6-t2.*t4.*t7.*t9.*t14.*t20.*1.27008e5+t3.*t5.*t7.*t9.*t12.*t20.*1.6128e4+t2.*t5.*t6.*t7.*t18.*t19.*1.6963392e7+t2.*t7.*t9.*t10.*t12.*t17.*2.8224e4+t3.*t5.*t7.*t9.*t16.*t17.*3.528e3+t4.*t6.*t7.*t9.*t11.*t20.*1.4112e4+t2.*t4.*t7.*t9.*t16.*t20.*4.2336e4-t2.*t6.*t9.*t10.*t11.*t20.*2.13696e5+t4.*t6.*t8.*t9.*t11.*t20.*4.18608e5+t2.*t6.*t8.*t11.*t12.*t20.*8.641728e6+t2.*t6.*t9.*t11.*t12.*t20.*5.58144e5+t6.*t7.*t8.*t9.*t11.*t20.*2.80476e5-t4.*t7.*t8.*t9.*t14.*t20.*9.5256e4+t5.*t6.*t7.*t8.*t18.*t19.*1.2722544e7+t7.*t8.*t9.*t10.*t12.*t17.*2.1168e4+t4.*t7.*t8.*t9.*t16.*t20.*3.1752e4-t6.*t8.*t9.*t10.*t11.*t20.*1.60272e5+t2.*t7.*t9.*t14.*t16.*t17.*4.2336e4+t6.*t7.*t9.*t11.*t12.*t20.*1.4112e4+t6.*t8.*t9.*t11.*t12.*t20.*4.18608e5+t3.*t6.*t7.*t14.*t18.*t19.*9.45189e5-t5.*t6.*t7.*t12.*t18.*t19.*3.60072e5-t6.*t9.*t10.*t11.*t14.*t17.*2.8224e4+t3.*t5.*t6.*t11.*t17.*t26.*3.360672e6-t6.*t9.*t10.*t11.*t12.*t20.*8.064e3-t3.*t6.*t7.*t16.*t18.*t19.*3.15063e5+t6.*t9.*t10.*t11.*t16.*t17.*9.408e3+t7.*t8.*t9.*t14.*t16.*t17.*3.1752e4+t6.*t9.*t11.*t14.*t16.*t20.*2.1168e4-t2.*t6.*t9.*t11.*t20.*t49.*1.0584e4-t2.*t6.*t9.*t11.*t20.*t50.*9.5256e4-t6.*t8.*t9.*t11.*t20.*t49.*7.938e3-t6.*t8.*t9.*t11.*t20.*t50.*7.1442e4+t2.*t3.*t5.*t6.*t8.*t11.*t17.*2.0164032e7+t2.*t3.*t5.*t6.*t9.*t11.*t17.*1.302336e6-t2.*t3.*t5.*t7.*t9.*t14.*t17.*3.1752e4+t2.*t3.*t5.*t7.*t9.*t12.*t20.*4.8384e4+t3.*t5.*t6.*t7.*t9.*t11.*t17.*3.2928e4+t2.*t3.*t5.*t7.*t9.*t16.*t17.*1.0584e4+t2.*t4.*t6.*t7.*t9.*t11.*t20.*4.2336e4+t3.*t5.*t6.*t8.*t9.*t11.*t17.*9.76752e5-t3.*t5.*t7.*t8.*t9.*t14.*t17.*2.3814e4+t3.*t5.*t7.*t8.*t9.*t12.*t20.*3.6288e4+t3.*t5.*t7.*t8.*t9.*t16.*t17.*7.938e3+t4.*t6.*t7.*t8.*t9.*t11.*t20.*3.1752e4+t2.*t6.*t7.*t9.*t11.*t12.*t20.*4.2336e4+t3.*t5.*t6.*t9.*t11.*t14.*t20.*9.072e3-t2.*t3.*t6.*t7.*t14.*t18.*t19.*5.041008e6+t2.*t5.*t6.*t7.*t12.*t18.*t19.*1.920384e6-t2.*t6.*t9.*t10.*t11.*t14.*t17.*8.4672e4-t2.*t6.*t9.*t10.*t11.*t12.*t20.*2.4192e4-t3.*t5.*t6.*t9.*t11.*t16.*t20.*3.024e3+t2.*t3.*t6.*t7.*t16.*t18.*t19.*1.680336e6+t2.*t6.*t9.*t10.*t11.*t16.*t17.*2.8224e4+t6.*t7.*t8.*t9.*t11.*t12.*t20.*3.1752e4-t3.*t6.*t7.*t8.*t14.*t18.*t19.*3.780756e6+t5.*t6.*t7.*t8.*t12.*t18.*t19.*1.440288e6-t6.*t8.*t9.*t10.*t11.*t14.*t17.*6.3504e4-t6.*t8.*t9.*t10.*t11.*t12.*t20.*1.8144e4+t3.*t6.*t7.*t8.*t16.*t18.*t19.*1.260252e6+t6.*t8.*t9.*t10.*t11.*t16.*t17.*2.1168e4+t2.*t6.*t9.*t11.*t14.*t16.*t20.*6.3504e4+t6.*t8.*t9.*t11.*t14.*t16.*t20.*4.7628e4+t2.*t3.*t5.*t6.*t7.*t9.*t11.*t17.*9.8784e4+t3.*t5.*t6.*t7.*t8.*t9.*t11.*t17.*7.4088e4+t2.*t3.*t5.*t6.*t9.*t11.*t14.*t20.*2.7216e4-t2.*t3.*t5.*t6.*t9.*t11.*t16.*t20.*9.072e3+t3.*t5.*t6.*t8.*t9.*t11.*t14.*t20.*2.0412e4-t3.*t5.*t6.*t8.*t9.*t11.*t16.*t20.*6.804e3).*(-1.0./5.0);t6.*t43.*t52.*t99.*1.0002e2];