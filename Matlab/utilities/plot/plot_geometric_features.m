function plot_geometric_features(p,fignum,features,lzid)
    % p - struct of parameters.
    % fignum - figure number
    % features - (1x13) cell which consists of features where each entry is one feature type.
    % lzid - laser id
    % p.eb_feats - true/false : enable/disable this figure of features
    % p.eb_clusteronly - true : plots features computed on only a cluster, false: plots features computed on entire scanline.

    [L,P,S,O,A,EigSum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,angle] = deal(features{:});
% [L,P,S,O,A,EigSum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA] = deal(features{:});
    N = size(L,1);
    pnt_index      = 1:N;    
    eb_feats       = p.eb_feats;  
    eb_clusteronly = p.eb_clusteronly;

    if eb_feats
        % subplot grid's row & col count.
        gr = 4; gc = 3;

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
    
        % Plot Omnivariance
        subplot(gr,gc,2); 
        hold on; grid on
        plot(pnt_index,O,'.')
        ylabel('Omnivariance')
        
        % Plot Anisotropy
        subplot(gr,gc,3);
        hold on; grid on
        plot(pnt_index,A,'.')
        xlabel('point index')
        title('Anisotropy')
        
        % % Plot Eigen Entropy
        % subplot(gr,gc,3); 
        % hold on; grid on
        % plot(pnt_index,Entropy,'.')
        % ylabel('EigEntropy')

        % Plot Surface Variation
        subplot(gr,gc,4); 
        hold on; grid on
        plot(pnt_index,SV,'.')
        ylabel('Surface Var.')

        % Plot Eigen value Sum
        subplot(gr,gc,5); 
        hold on; grid on
        plot(pnt_index,EigSum,'.')
        ylabel('EigenVal Sum')

        % Plot Smoothness
        subplot(gr,gc,6); hold on; grid on
        plot(pnt_index,Sm,'.')
        ylabel('Smoothness')
        
        % Plot Moments
        subplot(gr,gc,7);
        hold on; grid on
        plot(pnt_index,FOFA,'.')
        ylabel('1stOrd1stEVec')
    
        % Plot FOSA
        subplot(gr,gc,8)
        hold on; grid on
        plot(pnt_index,FOSA,'.')
        ylabel('1stOrd2ndEVec')
    
        % Plot SOFA
        subplot(gr,gc,9)
        hold on; grid on
        plot(pnt_index,SOFA,'.')
        xlabel('point index')
        ylabel('2ndOrd1stEVec')
    
        % Plot SOSA 
        subplot(gr,gc,10)
        hold on; grid on
        plot(pnt_index,SOSA,'.')
        xlabel('point index')
        ylabel('2ndOrd2ndEVec')

        % Plot Dot Product 
        subplot(gr,gc,11)
        hold on; grid on
        plot(pnt_index,dp,'.')
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

