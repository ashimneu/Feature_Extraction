function plot_geometric_features_v3(p,fignum,f,lzid)
    % p - struct of parameters.
    % fignum - figure number
    % f - struct of features
    % lzid - laser id
    % p.eb_feats - true/false : enable/disable this figure of features
    % p.eb_clusteronly - true : plots features computed on only a cluster, false: plots features computed on entire scanline.

    L = f.L;
    P = f.P;
    S = f.S;
    O = f.O;
    A = f.A;
    EigSum = f.Esum;
    SV = f.SV;
    Entropy = f.Entropy;
    Sm = f.Sm;
    FOFA = f.FOFA;
    FOSA = f.FOSA;
    SOFA = f.SOFA;
    SOSA = f.SOSA;
    dp = f.dotprod;
    alpha = f.angle;

    N = size(L,1);
    pnt_index      = 1:N;    
    eb_feats       = p.eb_feats;  
    eb_clusteronly = p.eb_clusteronly;

    if eb_feats
        % subplot grid's row & col count.
        gr=2; gc =2;

        % Figure of geometric features.
        figure(fignum); clf;

        % Plot Linearity        
        subplot(gr,gc,1); 
        hold on; grid on
        plot(pnt_index,L,'.')
        xlabel('point index')
        ylabel('Linearity')
        
        %  % Plot Planarity
        % subplot(gr,gc,2); 
        % hold on; grid on
        % plot(pnt_index,P,'.')
        % title('Planarity')
            
        % % Plot Sphericity
        % subplot(gr,gc,2); 
        % hold on; grid on
        % plot(pnt_index,S,'.')
        % title('Sphericity')
        
        % Plot angle (between principal directions 
        % of preceeding & succedding points)
        subplot(gr,gc,2); 
        hold on; grid on
        % flagged alpha (i.e. computation wasn't feasible)
        idx_flagged1 = find(alpha == -2);
        idx_flagged2 = find(alpha == -3);
        % Unflagged alpha (i.e. computation was feasible)
        idx_unflagged = find(alpha >= 0); 
        plot(pnt_index(idx_unflagged),alpha(idx_unflagged),'.')
        plot(pnt_index(idx_flagged1),alpha(idx_flagged1),'r.')
        plot(pnt_index(idx_flagged2),alpha(idx_flagged2),'g.')
        legend('nfl','fl','fl2','location','best')
        ylabel('Alpha')

        % % Plot Omnivariance
        % subplot(gr,gc,2); 
        % hold on; grid on
        % plot(pnt_index,O,'.')
        % ylabel('Omnivariance')
        
        % % Plot Anisotropy
        % subplot(gr,gc,3);
        % hold on; grid on
        % plot(pnt_index,A,'.')
        % xlabel('point index')
        % ylabel('Anisotropy')
        
        % % Plot Eigen Entropy
        % subplot(gr,gc,3); 
        % hold on; grid on
        % plot(pnt_index,Entropy,'.')
        % ylabel('EigEntropy')

        % % Plot Surface Variation
        % subplot(gr,gc,4); 
        % hold on; grid on
        % plot(pnt_index,SV,'.')
        % ylabel('Surface Var.')

        % Plot Eigen value Sum
        subplot(gr,gc,3); 
        hold on; grid on
        plot(pnt_index,EigSum,'.')
        ylabel('EigenVal Sum')

        % % Plot Smoothness
        % subplot(gr,gc,4); 
        % hold on; grid on
        % plot(pnt_index,Sm,'.')
        % ylabel('Smoothness')
        
        % % Plot Moments
        % subplot(gr,gc,7);
        % hold on; grid on
        % plot(pnt_index,FOFA,'.')
        % ylabel('1stOrd1stEVec')
        % 
        % % Plot FOSA
        % subplot(gr,gc,8)
        % hold on; grid on
        % plot(pnt_index,FOSA,'.')
        % ylabel('1stOrd2ndEVec')
        % 
        % % Plot SOFA
        % subplot(gr,gc,9)
        % hold on; grid on
        % plot(pnt_index,SOFA,'.')
        % xlabel('point index')
        % ylabel('2ndOrd1stEVec')
        % 
        % % Plot SOSA 
        % subplot(gr,gc,10)
        % hold on; grid on
        % plot(pnt_index,SOSA,'.')
        % xlabel('point index')
        % ylabel('2ndOrd2ndEVec')

        % Plot Dot Product        
        subplot(gr,gc,4)
        hold on; grid on
        % flagged dp (i.e. computation wasn't feasible)
        idx_flagged1 = find(dp == 2);
        idx_flagged2 = find(dp == 3);
        % Unflagged dp (i.e. computation was feasible)
        idx_unflagged = find(dp < 2);
        plot(pnt_index(idx_unflagged),dp(idx_unflagged),'.')
        plot(pnt_index(idx_flagged1),dp(idx_flagged1),'r.')
        plot(pnt_index(idx_flagged2),dp(idx_flagged2),'g.')
        legend('nfl','fl1','fl2','location','best')
        xlabel('point index')
        ylabel('dot product')
        
        if eb_clusteronly == true
            % title includes cluster number, sliding window size, # of
            % points in the cluster
            sgtitle({['Features computed on cluster #',num2str(p.clust_id),' only.'], ...
            ['window size=',num2str(p.sz_window), ', point count=',num2str(N)]});
        else
            % title only includes laser number.
            sgtitle({['Features, laser',num2str(lzid(1))]})
        end
    end
end

