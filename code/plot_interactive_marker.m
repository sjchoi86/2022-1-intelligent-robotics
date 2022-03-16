function plot_interactive_marker(varargin)

persistent h
global g_im
% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:10, h{i,j}.first_flag = true; end; end; end

% Parse input arguments
p = inputParser;
addParameter(p,'fig_idx',1);
addParameter(p,'subfig_idx',1);
addParameter(p,'T',eye(4,4));       % pose
addParameter(p,'clen',1.0);         % center length
addParameter(p,'clen_er',1.2);
addParameter(p,'tlw',4.0);          % translation line width
addParameter(p,'rlw',2.0);          % rotation line width
addParameter(p,'sr',0.1);           % sphere radius
addParameter(p,'fa',0.6);
addParameter(p,'USE_DRAGZOOM',true);
addParameter(p,'VERBOSE',false);
parse(p,varargin{:});
fig_idx         = p.Results.fig_idx;
subfig_idx      = p.Results.subfig_idx;
T               = p.Results.T;
clen            = p.Results.clen;
clen_er         = p.Results.clen_er;
tlw             = p.Results.tlw;
rlw             = p.Results.rlw;
sr              = p.Results.sr;
fa              = p.Results.fa;
USE_DRAGZOOM    = p.Results.USE_DRAGZOOM;
VERBOSE         = p.Results.VERBOSE;

% Helper functions for interactive marker
% ---------------------- X-translation ----------------------
    function bdf_xsphere(hcbo,evt) % x-translation click
        if VERBOSE
            fprintf('[%d-%d] x-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_xsphere;
        fig.WindowButtonUpFcn = @wbuf_xsphere;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[1,0,0]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_xsphere(hcbo,evt) % x-translation move
        if VERBOSE
            fprintf('[%d-%d] x-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 - clen_er * clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_xsphere(hcbo,evt) % x-translation unclick
        if VERBOSE
            fprintf('[%d-%d] x-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end
    function bdf_xsphere2(hcbo,evt) % x-translation click
        if VERBOSE
            fprintf('[%d-%d] x-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_xsphere2;
        fig.WindowButtonUpFcn = @wbuf_xsphere2;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[1,0,0]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_xsphere2(hcbo,evt) % x-translation move
        if VERBOSE
            fprintf('[%d-%d] x-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 + clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_xsphere2(hcbo,evt) % x-translation unclick
        if VERBOSE
            fprintf('[%d-%d] x-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end

% ---------------------- Y-translation ----------------------
    function bdf_ysphere(hcbo,evt) % y-translation click
        if VERBOSE
            fprintf('[%d-%d] y-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_ysphere;
        fig.WindowButtonUpFcn = @wbuf_ysphere;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[0,1,0]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_ysphere(hcbo,evt) % y-translation move
        if VERBOSE
            fprintf('[%d-%d] y-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 - clen_er * clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_ysphere(hcbo,evt) % y-translation unclick
        if VERBOSE
            fprintf('[%d-%d] y-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end
    function bdf_ysphere2(hcbo,evt) % y-translation click
        if VERBOSE
            fprintf('[%d-%d] y-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_ysphere2;
        fig.WindowButtonUpFcn = @wbuf_ysphere2;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[0,1,0]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_ysphere2(hcbo,evt) % y-translation move
        if VERBOSE
            fprintf('[%d-%d] y-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 + clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_ysphere2(hcbo,evt) % y-translation unclick
        if VERBOSE
            fprintf('[%d-%d] y-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end

% ---------------------- Z-translation ----------------------
    function bdf_zsphere(hcbo,evt) % z-translation click
        if VERBOSE
            fprintf('[%d-%d] z-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_zsphere;
        fig.WindowButtonUpFcn = @wbuf_zsphere;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[0,0,1]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_zsphere(hcbo,evt) % z-translation move
        if VERBOSE
            fprintf('[%d-%d] z-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 - clen_er * clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_zsphere(hcbo,evt) % z-translation unclick
        if VERBOSE
            fprintf('[%d-%d] z-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end
    function bdf_zsphere2(hcbo,evt) % z-translation click
        if VERBOSE
            fprintf('[%d-%d] z-translation click. \n',fig_idx,subfig_idx);
        end
        if 0
            h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        else
            h{fig_idx,subfig_idx}.ax = hcbo.Parent.Parent;
        end
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_zsphere2;
        fig.WindowButtonUpFcn = @wbuf_zsphere2;
        % Get direction?
        dir_temp = h{fig_idx,subfig_idx}.T(1:3,1:3)*[0,0,1]';
        h{fig_idx,subfig_idx}.dir2 = dir_temp / norm(dir_temp);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_zsphere2(hcbo,evt) % z-translation move
        if VERBOSE
            fprintf('[%d-%d] z-translation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        p1 = cp(1,:)';
        p2 = cp(2,:)';
        dir1 = p1 - p2;
        dir2 = h{fig_idx,subfig_idx}.dir2;
        n = cr3(dir1)*dir2;
        n1 = cr3(dir1)*n;
        origin2 = h{fig_idx,subfig_idx}.T(1:3,4);
        c2 = origin2 + ( ((p1 - origin2)'*n1 )/(dir2'*n1)) * dir2;
        p = c2 + clen * dir2;
        h{fig_idx,subfig_idx}.T(1:3,4) = p;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_zsphere2(hcbo,evt) % z-translation unclick
        if VERBOSE
            fprintf('[%d-%d] z-translation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end

% ---------------------- X-rotation ----------------------
    function bdf_xrot(hcbo,evt) % x-rotation button click
        if VERBOSE
            fprintf('[%d-%d] x-rotation click. \n',fig_idx,subfig_idx);
        end
        x = hcbo.XData;
        y = hcbo.YData;
        z = hcbo.ZData;
        P = [x;y;z];
        p1 = P(:,1);
        p2 = P(:,2);
        p3 = P(:,3);
        sn_temp = cr3(p1-p2)*(p2-p3);
        h{fig_idx,subfig_idx}.sn = sn_temp/norm(sn_temp);
        h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        h{fig_idx,subfig_idx}.apoc = p1;
        s = nume/deno;
        h{fig_idx,subfig_idx}.intersecini = q1 + s*(q2-q1);
        qa = P(:,1);
        qb = P(:,33);
        h{fig_idx,subfig_idx}.cc = 0.5*(qa +qb);
        h{fig_idx,subfig_idx}.R0 = h{fig_idx,subfig_idx}.T(1:3, 1:3);
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_xrot;
        fig.WindowButtonUpFcn = @wbuf_xrot;
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_xrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] x-rotation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        p1 = h{fig_idx,subfig_idx}.apoc;
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        s = nume/deno;
        current_intersec = q1 + s*(q2-q1);
        aa = (current_intersec - h{fig_idx,subfig_idx}.cc);
        aa = aa/norm(aa);
        bb = (h{fig_idx,subfig_idx}.intersecini - h{fig_idx,subfig_idx}.cc);
        bb =  bb/norm(bb);
        costheta = aa'*bb;
        st = -cr3(aa)*bb;
        for i_idx = 1:3
            if abs(h{fig_idx,subfig_idx}.sn(i_idx)) > 1e-10
                sintheta = st(i_idx)/h{fig_idx,subfig_idx}.sn(i_idx);
                break
            end
        end
        theta = atan2(sintheta, costheta);
        R = expm(theta*cr3(h{fig_idx,subfig_idx}.sn));
        h{fig_idx,subfig_idx}.T(1:3, 1:3) = R*h{fig_idx,subfig_idx}.R0;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_xrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] x-rotation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end

% ---------------------- Y-rotation ----------------------
    function bdf_yrot(hcbo,evt) % y-rotation button click
        if VERBOSE
            fprintf('[%d-%d] y-rotation click. \n',fig_idx,subfig_idx);
        end
        x = hcbo.XData;
        y = hcbo.YData;
        z = hcbo.ZData;
        P = [x;y;z];
        p1 = P(:,1);
        p2 = P(:,2);
        p3 = P(:,3);
        sn_temp = cr3(p1-p2)*(p2-p3);
        h{fig_idx,subfig_idx}.sn = sn_temp/norm(sn_temp);
        h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        h{fig_idx,subfig_idx}.apoc = p1;
        s = nume/deno;
        h{fig_idx,subfig_idx}.intersecini = q1 + s*(q2-q1);
        qa = P(:,1);
        qb = P(:,33);
        h{fig_idx,subfig_idx}.cc = 0.5*(qa + qb);
        h{fig_idx,subfig_idx}.R0 = h{fig_idx,subfig_idx}.T(1:3, 1:3);
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_yrot;
        fig.WindowButtonUpFcn = @wbuf_yrot;
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_yrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] y-rotation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        p1 = h{fig_idx,subfig_idx}.apoc;
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        s = nume/deno;
        current_intersec = q1 + s*(q2-q1);
        aa = (current_intersec - h{fig_idx,subfig_idx}.cc);
        aa = aa/norm(aa);
        bb = (h{fig_idx,subfig_idx}.intersecini - h{fig_idx,subfig_idx}.cc);
        bb =  bb/norm(bb);
        costheta = aa'*bb;
        st = -cr3(aa)*bb;
        for i_idx = 1:3
            if abs(h{fig_idx,subfig_idx}.sn(i_idx)) > 1e-10
                sintheta = st(i_idx)/h{fig_idx,subfig_idx}.sn(i_idx);
                break
            end
        end
        theta = atan2(sintheta, costheta);
        R = expm(theta*cr3(h{fig_idx,subfig_idx}.sn));
        h{fig_idx,subfig_idx}.T(1:3, 1:3) = R*h{fig_idx,subfig_idx}.R0;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_yrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] y-rotation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end

% ---------------------- Z-rotation ----------------------
    function bdf_zrot(hcbo,evt) % z-rotation button click
        if VERBOSE
            fprintf('[%d-%d] z-rotation click. \n',fig_idx,subfig_idx);
        end
        x = hcbo.XData;
        y = hcbo.YData;
        z = hcbo.ZData;
        P = [x;y;z];
        p1 = P(:,1);
        p2 = P(:,2);
        p3 = P(:,3);
        sn_temp = cr3(p1-p2)*(p2-p3);
        h{fig_idx,subfig_idx}.sn = sn_temp/norm(sn_temp);
        h{fig_idx,subfig_idx}.ax = hcbo.Parent;
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        h{fig_idx,subfig_idx}.apoc = p1;
        s = nume/deno;
        h{fig_idx,subfig_idx}.intersecini = q1 + s*(q2-q1);
        qa = P(:,1);
        qb = P(:,33);
        h{fig_idx,subfig_idx}.cc = 0.5*(qa + qb);
        h{fig_idx,subfig_idx}.R0 = h{fig_idx,subfig_idx}.T(1:3, 1:3);
        fig = h{fig_idx,subfig_idx}.ax.Parent;
        fig.WindowButtonMotionFcn = @wbmf_zrot;
        fig.WindowButtonUpFcn = @wbuf_zrot;
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'clicked';
    end
    function wbmf_zrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] z-rotation move. \n',fig_idx,subfig_idx);
        end
        cp = h{fig_idx,subfig_idx}.ax.CurrentPoint;
        q1 = cp(1,:)';
        q2 = cp(2,:)';
        p1 = h{fig_idx,subfig_idx}.apoc;
        nume = -h{fig_idx,subfig_idx}.sn'*(q1 - p1);
        deno = h{fig_idx,subfig_idx}.sn'*(q2 - q1);
        s = nume/deno;
        current_intersec = q1 + s*(q2-q1);
        aa = (current_intersec - h{fig_idx,subfig_idx}.cc);
        aa = aa/norm(aa);
        bb = (h{fig_idx,subfig_idx}.intersecini - h{fig_idx,subfig_idx}.cc);
        bb =  bb/norm(bb);
        costheta = aa'*bb;
        st = -cr3(aa)*bb;
        for i_idx = 1:3
            if abs(h{fig_idx,subfig_idx}.sn(i_idx)) > 1e-10
                sintheta = st(i_idx)/h{fig_idx,subfig_idx}.sn(i_idx);
                break
            end
        end
        theta = atan2(sintheta, costheta);
        R = expm(theta*cr3(h{fig_idx,subfig_idx}.sn));
        h{fig_idx,subfig_idx}.T(1:3, 1:3) = R*h{fig_idx,subfig_idx}.R0;
        % Update plot
        update_plot(h{fig_idx,subfig_idx}.T);
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updating';
    end
    function wbuf_zrot(hcbo,evt)
        if VERBOSE
            fprintf('[%d-%d] z-rotation unclick. \n',fig_idx,subfig_idx);
        end
        % Clear button motion and button up function handlers
        set(hcbo,'WindowButtonMotionFcn','')
        set(hcbo,'WindowButtonUpFcn','')
        if USE_DRAGZOOM
            dragzoom; % go back to use dragzoom
        end
        % Set status
        g_im{fig_idx,subfig_idx}.status = 'updated';
    end


    function pcross = cr3( p )
        pcross = [0, -p(3), p(2); p(3), 0, -p(1); -p(2), p(1), 0];
    end

    function [Xs,Ys,Zs,X_xtrn,X_xtrn2,X_ytrn,X_ytrn2,X_ztrn,X_ztrn2,X_xrot,X_yrot,X_zrot,...
            xyz_fv,xyz_fv_small] = ...
            get_ingredients(T)
        p = t2p(T)';
        R = t2r(T);
        % Plot interactive marker
        theta = (0:pi/32:2*pi)';
        x = clen*cos(theta);
        y = clen*sin(theta);
        l = length(x);
        z = zeros(l,1);
        if 0
            [Xs, Ys, Zs] = sphere(40);
        else
            [Xs, Ys, Zs] = ellipsoid(0,0,0,1,1,1,20);
            xyz_fv = surf2patch(sr*Xs, sr*Ys, sr*Zs);
            xyz_fv_small = surf2patch(sr*Xs*0.6,sr*Ys*0.6,sr*Zs*0.6);
        end
        Xs = sr*Xs;
        Ys = sr*Ys;
        Zs = sr*Zs;
        X_xtrn = [-clen, 0, 0; clen_er*clen, 0, 0];
        X_xtrn = X_xtrn*R' + ones(2,1)*p;
        X_xtrn2 = [-0, 0, 0; clen_er*clen, 0, 0];
        X_xtrn2 = X_xtrn2*R' + ones(2,1)*p;
        X_ytrn = [0, -clen, 0; 0, clen_er*clen, 0];
        X_ytrn = X_ytrn*R' + ones(2,1)*p;
        X_ytrn2 = [0, -0, 0; 0, clen_er*clen, 0];
        X_ytrn2 = X_ytrn2*R' + ones(2,1)*p;
        X_ztrn = [0, 0, -clen; 0, 0, clen_er*clen];
        X_ztrn = X_ztrn*R' + ones(2,1)*p;
        X_ztrn2 = [0, 0, -0; 0, 0, clen_er*clen];
        X_ztrn2 = X_ztrn2*R' + ones(2,1)*p;
        X_xrot = [z,x,y]*R' + ones(l,1)*p;
        X_yrot = [y,z,x]*R' + ones(l,1)*p;
        X_zrot = [x,y,z]*R' + ones(l,1)*p;
    end

    function update_plot(T)
        [Xs,Ys,Zs,X_xtrn,X_xtrn2,X_ytrn,X_ytrn2,X_ztrn,X_ztrn2,X_xrot,X_yrot,X_zrot,...
            xyz_fv,xyz_fv_small] = ...
            get_ingredients(T);
        % x-translation
        h{fig_idx,subfig_idx}.xaxis.XData = X_xtrn2(:,1);
        h{fig_idx,subfig_idx}.xaxis.YData = X_xtrn2(:,2);
        h{fig_idx,subfig_idx}.xaxis.ZData = X_xtrn2(:,3);
        if 0
            h{fig_idx,subfig_idx}.xsphere.XData = Xs+X_xtrn(2,1);
            h{fig_idx,subfig_idx}.xsphere.YData = Ys+X_xtrn(2,2);
            h{fig_idx,subfig_idx}.xsphere.ZData = Zs+X_xtrn(2,3);
            h{fig_idx,subfig_idx}.xsphere2.XData = 0.5*Xs+X_xtrn(1,1);
            h{fig_idx,subfig_idx}.xsphere2.YData = 0.5*Ys+X_xtrn(1,2);
            h{fig_idx,subfig_idx}.xsphere2.ZData = 0.5*Zs+X_xtrn(1,3);
        else
            tform = p2t(X_xtrn(2,:));
            set(h{fig_idx,subfig_idx}.xsphere_t,'Matrix',tform);
            tform = p2t(X_xtrn(1,:));
            set(h{fig_idx,subfig_idx}.xsphere2_t,'Matrix',tform);
        end
        % y-translation
        h{fig_idx,subfig_idx}.yaxis.XData = X_ytrn2(:,1);
        h{fig_idx,subfig_idx}.yaxis.YData = X_ytrn2(:,2);
        h{fig_idx,subfig_idx}.yaxis.ZData = X_ytrn2(:,3);
        if 0
            h{fig_idx,subfig_idx}.ysphere.XData = Xs+X_ytrn(2,1);
            h{fig_idx,subfig_idx}.ysphere.YData = Ys+X_ytrn(2,2);
            h{fig_idx,subfig_idx}.ysphere.ZData = Zs+X_ytrn(2,3);
            h{fig_idx,subfig_idx}.ysphere2.XData = 0.5*Xs+X_ytrn(1,1);
            h{fig_idx,subfig_idx}.ysphere2.YData = 0.5*Ys+X_ytrn(1,2);
            h{fig_idx,subfig_idx}.ysphere2.ZData = 0.5*Zs+X_ytrn(1,3);
        else
            tform = p2t(X_ytrn(2,:));
            set(h{fig_idx,subfig_idx}.ysphere_t,'Matrix',tform);
            tform = p2t(X_ytrn(1,:));
            set(h{fig_idx,subfig_idx}.ysphere2_t,'Matrix',tform);
        end
        % z-translation
        h{fig_idx,subfig_idx}.zaxis.XData = X_ztrn2(:,1);
        h{fig_idx,subfig_idx}.zaxis.YData = X_ztrn2(:,2);
        h{fig_idx,subfig_idx}.zaxis.ZData = X_ztrn2(:,3);
        if 0
            h{fig_idx,subfig_idx}.zsphere.XData = Xs+X_ztrn(2,1);
            h{fig_idx,subfig_idx}.zsphere.YData = Ys+X_ztrn(2,2);
            h{fig_idx,subfig_idx}.zsphere.ZData = Zs+X_ztrn(2,3);
            h{fig_idx,subfig_idx}.zsphere2.XData = 0.5*Xs+X_ztrn(1,1);
            h{fig_idx,subfig_idx}.zsphere2.YData = 0.5*Ys+X_ztrn(1,2);
            h{fig_idx,subfig_idx}.zsphere2.ZData = 0.5*Zs+X_ztrn(1,3);
        else
            tform = p2t(X_ztrn(2,:));
            set(h{fig_idx,subfig_idx}.zsphere_t,'Matrix',tform);
            tform = p2t(X_ztrn(1,:));
            set(h{fig_idx,subfig_idx}.zsphere2_t,'Matrix',tform);
        end
        % x-rotation
        h{fig_idx,subfig_idx}.xrot.XData = X_xrot(:,1);
        h{fig_idx,subfig_idx}.xrot.YData = X_xrot(:,2);
        h{fig_idx,subfig_idx}.xrot.ZData = X_xrot(:,3);
        % y-rotation
        h{fig_idx,subfig_idx}.yrot.XData = X_yrot(:,1);
        h{fig_idx,subfig_idx}.yrot.YData = X_yrot(:,2);
        h{fig_idx,subfig_idx}.yrot.ZData = X_yrot(:,3);
        % z-rotation
        h{fig_idx,subfig_idx}.zrot.XData = X_zrot(:,1);
        h{fig_idx,subfig_idx}.zrot.YData = X_zrot(:,2);
        h{fig_idx,subfig_idx}.zrot.ZData = X_zrot(:,3);
        drawnow limitrate;
        % Update global interactive marker T
        g_im{fig_idx,subfig_idx}.T = h{fig_idx,subfig_idx}.T;
    end

if h{fig_idx,subfig_idx}.first_flag || ... % first flag
        ~ishandle(h{fig_idx,subfig_idx}.fig)
    h{fig_idx,subfig_idx}.first_flag = false;
    % Initialize figure
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    % Initialize T
    h{fig_idx,subfig_idx}.T = T;
    g_im{fig_idx,subfig_idx}.T = T;
    [Xs,Ys,Zs,X_xtrn,X_xtrn2,X_ytrn,X_ytrn2,X_ztrn,X_ztrn2,X_xrot,X_yrot,X_zrot,...
        xyz_fv,xyz_fv_small] = ...
        get_ingredients(T);
    
    % x-translation (axis and sphere)
    h{fig_idx,subfig_idx}.xaxis = line(X_xtrn2(:,1), X_xtrn2(:,2), X_xtrn2(:,3),...
        'color','r', 'LineWidth',tlw,'ButtonDownFcn','');
    if 0
        h{fig_idx,subfig_idx}.xsphere = surface(Xs+X_xtrn(2,1), Ys+X_xtrn(2,2), Zs+X_xtrn(2,3),...
            'facecolor','r', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_xsphere);
        h{fig_idx,subfig_idx}.xsphere2 = surface(0.5*Xs+X_xtrn(1,1), 0.5*Ys+X_xtrn(1,2), 0.5*Zs+X_xtrn(1,3),...
            'facecolor','r', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_xsphere2);
    else
        h{fig_idx,subfig_idx}.xsphere = patch(xyz_fv,...
            'FaceColor','r','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_xsphere);
        h{fig_idx,subfig_idx}.xsphere_t = hgtransform;
        set(h{fig_idx,subfig_idx}.xsphere,'parent',h{fig_idx,subfig_idx}.xsphere_t);
        tform = p2t(X_xtrn(2,:));
        set(h{fig_idx,subfig_idx}.xsphere_t,'Matrix',tform);
        h{fig_idx,subfig_idx}.xsphere2 = patch(xyz_fv_small,...
            'FaceColor','r','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_xsphere2);
        h{fig_idx,subfig_idx}.xsphere2_t = hgtransform;
        set(h{fig_idx,subfig_idx}.xsphere2,'parent',h{fig_idx,subfig_idx}.xsphere2_t);
        tform = p2t(X_xtrn(1,:));
        set(h{fig_idx,subfig_idx}.xsphere2_t,'Matrix',tform);
    end
    
    % y-translation
    h{fig_idx,subfig_idx}.yaxis = line(X_ytrn2(:,1), X_ytrn2(:,2), X_ytrn2(:,3),...
        'color','g', 'LineWidth',tlw,'ButtonDownFcn','');
    if 0
        h{fig_idx,subfig_idx}.ysphere = surface(Xs+X_ytrn(2,1), Ys+X_ytrn(2,2), Zs+X_ytrn(2,3),...
            'facecolor','g', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_ysphere);
        h{fig_idx,subfig_idx}.ysphere2 = surface(0.5*Xs+X_ytrn(1,1), 0.5*Ys+X_ytrn(1,2), 0.5*Zs+X_ytrn(1,3),...
            'facecolor','g', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_ysphere2);
    else
        h{fig_idx,subfig_idx}.ysphere = patch(xyz_fv,...
            'FaceColor','g','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_ysphere);
        h{fig_idx,subfig_idx}.ysphere_t = hgtransform;
        set(h{fig_idx,subfig_idx}.ysphere,'parent',h{fig_idx,subfig_idx}.ysphere_t);
        tform = p2t(X_ytrn(2,:));
        set(h{fig_idx,subfig_idx}.ysphere_t,'Matrix',tform);
        h{fig_idx,subfig_idx}.ysphere2 = patch(xyz_fv_small,...
            'FaceColor','g','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_ysphere2);
        h{fig_idx,subfig_idx}.ysphere2_t = hgtransform;
        set(h{fig_idx,subfig_idx}.ysphere2,'parent',h{fig_idx,subfig_idx}.ysphere2_t);
        tform = p2t(X_ytrn(1,:));
        set(h{fig_idx,subfig_idx}.ysphere2_t,'Matrix',tform);
    end
    
    % z-translation
    h{fig_idx,subfig_idx}.zaxis = line(X_ztrn2(:,1), X_ztrn2(:,2), X_ztrn2(:,3),...
        'color','b', 'LineWidth',tlw,'ButtonDownFcn','');
    if 0
        h{fig_idx,subfig_idx}.zsphere = surface(Xs+X_ztrn(2,1), Ys+X_ztrn(2,2), Zs+X_ztrn(2,3),...
            'facecolor','b', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_zsphere);
        h{fig_idx,subfig_idx}.zsphere2 = surface(0.5*Xs+X_ztrn(1,1), 0.5*Ys+X_ztrn(1,2), 0.5*Zs+X_ztrn(1,3),...
            'facecolor','b', 'linestyle','none','facelighting','gouraud','ButtonDownFcn',@bdf_zsphere2);
    else
        h{fig_idx,subfig_idx}.zsphere = patch(xyz_fv,...
            'FaceColor','b','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_zsphere);
        h{fig_idx,subfig_idx}.zsphere_t = hgtransform;
        set(h{fig_idx,subfig_idx}.zsphere,'parent',h{fig_idx,subfig_idx}.zsphere_t);
        tform = p2t(X_ztrn(2,:));
        set(h{fig_idx,subfig_idx}.zsphere_t,'Matrix',tform);
        h{fig_idx,subfig_idx}.zsphere2 = patch(xyz_fv_small,...
            'FaceColor','b','FaceAlpha',fa,'EdgeColor','none',...
            'facelighting','gouraud','ButtonDownFcn',@bdf_zsphere2);
        h{fig_idx,subfig_idx}.zsphere2_t = hgtransform;
        set(h{fig_idx,subfig_idx}.zsphere2,'parent',h{fig_idx,subfig_idx}.zsphere2_t);
        tform = p2t(X_ztrn(1,:));
        set(h{fig_idx,subfig_idx}.zsphere2_t,'Matrix',tform);
    end
    
    % x-rotation
    h{fig_idx,subfig_idx}.xrot = line(X_xrot(:,1), X_xrot(:,2), X_xrot(:,3),...
        'color','r', 'LineWidth',rlw,'ButtonDownFcn',@bdf_xrot);
    
    % y-rotation
    h{fig_idx,subfig_idx}.yrot = line(X_yrot(:,1), X_yrot(:,2), X_yrot(:,3),...
        'color','g', 'LineWidth',rlw,'ButtonDownFcn',@bdf_yrot);
    
    % z-rotation
    h{fig_idx,subfig_idx}.zrot = line(X_zrot(:,1), X_zrot(:,2), X_zrot(:,3),...
        'color','b', 'LineWidth',rlw,'ButtonDownFcn',@bdf_yrot);
    
    % Set status
    g_im{fig_idx,subfig_idx}.status = 'initialized';
    
    
else
    % Update T
    if VERBOSE
        fprintf('[%d-%d] update T. \n',fig_idx,subfig_idx);
    end
    h{fig_idx,subfig_idx}.T = T;
    update_plot(h{fig_idx,subfig_idx}.T);
    
    % Set status
    g_im{fig_idx,subfig_idx}.status = 'updated';
end
end