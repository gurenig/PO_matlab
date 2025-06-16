%> @file mycross.m
%> @brief Computes the cross product of two 3D vector fields component-wise.
%>
%> Given the Cartesian components of two 3D vectors or vector fields A and B,
%> this function computes their cross product C = A x B in Cartesian form.
%>
%> @param Ax X-component of vector A
%> @param Ay Y-component of vector A
%> @param Az Z-component of vector A
%> @param Bx X-component of vector B
%> @param By Y-component of vector B
%> @param Bz Z-component of vector B
%>
%> @retval Cx X-component of the cross product A x B
%> @retval Cy Y-component of the cross product A x B
%> @retval Cz Z-component of the cross product A x B
function [Cx,Cy,Cz] = mycross(Ax,Ay,Az,Bx,By,Bz)
    Cx = Ay.*Bz - Az.*By;
    Cy = Az.*Bx - Ax.*Bz;
    Cz = Ax.*By - Ay.*Bx;
end
