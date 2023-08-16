function pts = test_grow_random_poly(n,max_width,max_height)
%TEST_GROW_RANDOM_POLY Builds a random (concave) polygon.
% Uses a crystallization process to grow a polygon.
% takes optional arguments for number of sides (default = 10) and
% maximum allowable width and height (default=1)



if ~exist("n","var")
    n = 10;
end
if ~exist("max_width","var")
    max_width = 1;
end
if ~exist("max_height","var")
    max_height = 1;
end

%nucleation triangle
pts = zeros(n,2);
pts(1:3,1) = rand(3,1) * max_width;
pts(1:3,2) = rand(3,1) * max_height;

%ensure we have CCW motion: 2d cross product
u = pts(2,:) - pts(1,:);
v = pts(3,:) - pts(1,:);
if u(1)*v(2) - u(2)*v(1) < 0
    %is CW, flip last two pts
    pts(2:3,:) = pts([3,2],:);
end


%max attempts before trying a new edge
max_attempts_edge = 1000;
%max attempts before defaulting to the middle of the last tried edge
max_attempts_pt = 10;

for k=4:n
    %choose point
    found_pt = 0;
    for pt_attmpt=1:max_attempts_pt
        edge = randi(k-1);
        a = pts(edge-1 + (k-1)*(edge == 1),:);
        b = pts(edge,:);
        normal = flip(b-a) .* [1,-1];
        normal = normal ./ vecnorm(normal);
        
        %for this edge; try points
        for e_attmpt=1:max_attempts_edge
            p = zeros(1,2);
            p(1) = rand() * max_width;
            p(2) = rand() * max_height;
            
            if (dot(normal, p-a) <= 0) ||...
                    intersects(pts(1:(k-1),:),a,p) ||...
                    intersects(pts(1:(k-1),:),p,b)
                continue
            end
            %valid point
            found_pt = 1;
            break;

        end

        %did we find a good point? if no; keep going
        if found_pt
            break
        end
    end
    
    if ~found_pt
        %no point found; just choose a point on this edge
        p = a + (b-a)*(0.5 + 0.75*(rand() - 0.5));
    end

    %push p in between a and b
    pts((edge+1):k,:) = pts(edge:(k-1),:);
    pts(edge,:) = p;

%     for edge=1:k
%         a = pts(edge-1 + k*(edge == 1),:);
%         b = pts(edge,:);
%         if intersects(pts(1:k,:),a,b)
%             keyboard
%         end
%     end
end

end

function [does] = intersects(poly,a,b)
% returns 1 if the interior of ab intersects the polygon (boundary), or 0
% otherwise.
%Note: for our purposes, since a or b will be on the boundary,
%this check is sufficient to tell if ab intersects the interior of poly.

EPS = 1e-10;

u = b-a;

%segments: p_k - p_{k-1}
v = circshift(poly,1,1) - poly;

%calc intersect point
for i=1:size(poly,1)
    %a+su = p+tv
    collis = [-u;v(i,:)]' \ (a-poly(i,:))';
    
    %interior of segment a->b; closure of p->p+v
    if EPS <= collis(1) && collis(1) <= 1-EPS &&...
         0 <= collis(2) && collis(2) <= 1

        does = 1;
        return
    end

end


does = 0;
end