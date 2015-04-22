% scatter 3D of the manifold for a 2x2 covariance matrices set
%
function manifold_plot(COV,label,boundary)
    Nt = size(COV,3);
    if nargin<2, label =  ones(Nt,1); end
    if nargin<3, boundary = 0; end
    
    class_label = unique(label);
    N_label = length(class_label);
    clrs = hsv(N_label);
    
    color = zeros(Nt,3);
    for i=1:N_label
        ix = label==class_label(i);
        color(ix,:) = repmat(clrs(i,:),sum(ix),1);
    end
    x = squeeze(COV(1,1,:));
    y = squeeze(COV(2,2,:));
    z = squeeze(COV(1,2,:));
    scatter3(x,y,z,[],color,'filed','MarkerEdgeColor','k');
    grid on;
    view(60,25);

    if boundary
        mx = max(x);
        my = max(y);
        a = [0:mx/20:mx];
        c = [0:my/20:my];
        b = sqrt(a'*c);
        hold on;
        mesh(a,c,b','FaceColor','w','EdgeColor','k');
        mesh(a,c,-b','FaceColor','w','EdgeColor','k');hidden off;
    end