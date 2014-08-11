function Yd = getYd(Q)
% Input scheme:
%     Q(1:2) = qd; Q(3:4) = qd_dot; Q(5:6) = qd_dot_dot;
q1 = Q(1);
q2 = Q(2);
q_1 = Q(3);
q_2 = Q(4);
q__1 = Q(5);
q__2 = Q(6);

Yd1col = [q__1; 0];
Yd2col = [q__2; q__1+q__2];
Yd3col = [2*cos(q2)*q__1 + cos(q2)*q__2 - sin(q2)*q_1*q_2 - sin(q2)*(q_1+q_2)*q_2;...
  cos(q2)*q__1 + sin(q2)*q_1^2];
Yd4col = [q_1; 0];
Yd5col = [0;q_2];

Yd = [Yd1col, Yd2col, Yd3col, Yd4col, Yd5col];