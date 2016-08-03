function [T] = form_random_tetrads(X, Y, Z, dom);

T = -ones(5,1);

% Max Shape Factor
sfMax = 0.16025;

% Superpose periodicity to make search easier
% TODO: explicitly check periodicity
N = dom.N;
X = [X; X + dom.xl; X - dom.xl];
Y = [Y; Y + dom.yl; Y - dom.yl];
Z = [Z; Z + dom.zl; Z - dom.xl];
%Z = [Z; Z         ; Z         ];
Nper = 3*N;

% Randomize
np1 = unique(randi(dom.N, 1, ceil(dom.N/2)));
np2 = unique(randi(dom.N, 1, ceil(Nper/2)));
np3 = unique(randi(dom.N, 1, ceil(Nper/2)));
np4 = unique(randi(dom.N, 1, ceil(Nper/2)));

% Loop over ~all particles and find EW,NS,TB coordinates
for n1 = np1;
  % Set target particle
  p1.n = n1;
  p1.X = [X(n1); Y(n1); Z(n1)];

  for n2 = np2
    % make sure its not the 'same' particle, flipped over
    nmod = mod(n2,N);
    if nmod == 0
      nmod = N;
    end
    if (nmod == n1)
      continue;
    end

    % Pull target part info
    p2.n = nmod;
    p2.X = [X(n2); Y(n2); Z(n2)];

    for n3 = np3
      % make sure its not the 'same' particle, flipped over
      nmod = mod(n3,N);
      if nmod == 0
        nmod = N;
      end
      if (nmod == n1 | nmod == n2)
        continue;
      end

      % Pull target part info
      p3.n = nmod;
      p3.X = [X(n3); Y(n3); Z(n3)];

      for n4 = np4
        % make sure its not the 'same' particle, flipped over
        nmod = mod(n4,N);
        if nmod == 0
          nmod = N;
        end
        if (nmod == n1 | nmod == n2 | nmod == n3)
          continue;
        end

        % Pull target part info
        p4.n = nmod;
        p4.X = [X(n4); Y(n4); Z(n4)];

        %% Calculate initial shape factor
        rho_1 = (p2.X - p1.X)./sqrt(2);
        rho_2 = (2*p3.X - p2.X - p1.X)./sqrt(6);
        rho_3 = (3*p4.X - p3.X - p2.X - p1.X)./sqrt(12);

        g = [rho_1, rho_2, rho_3];
        G = g*transpose(g);
        lambda = eig(G);
        shapeFactor = (prod(lambda))^(1/3)/3/sum(lambda)/sfMax;

        %if shapeFactor > 2/3
          temp = [p1.n; p2.n; p3.n; p4.n; shapeFactor];
          T = [T, temp];
        %end
      end
    end
  end
end

% Remove initialization
if size(T,2) > 1
  T = T(:,2:end);
end
% sort columns of T from small to large
T = sort(T, 1);
% find only unique entries
% array gymnastics required due to no 'columns' option
T = unique(T', 'rows')';
