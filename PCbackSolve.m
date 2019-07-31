function y = PCbackSolve(x,A,L1,U1,P,Q,L2,U2,a1)
        x((a1+1):end) = (U2\(L2\x((a1+1):end)));
        x(1:a1) = x(1:a1) - A*x((a1+1):end);
        x(1:a1) = Q*(U1\(L1\(P*x(1:a1))));
        y = x;
end