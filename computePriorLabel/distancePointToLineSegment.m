

function distance = distancePointToLineSegment(p, a, b)
    % p is a point
    % a is the start of a line segment
    % b is the end of a line segment

    pa_distance = sum((a-p).^2);
    pb_distance = sum((b-p).^2);

    unitab = (a-b) / norm(a-b);

    if pa_distance < pb_distance
        d = dot((p-a)/norm(p-a), -unitab);

        if d > 0
            distance = sqrt(pa_distance * (1 - d*d));
        else
            distance = sqrt(pa_distance);
        end
    else
        d = dot((p-b)/norm(p-b), unitab);

        if d > 0
            distance = sqrt(pb_distance * (1 - d*d));
        else
            distance = sqrt(pb_distance);
        end
    end
end