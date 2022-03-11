% close all;

spectfile = 't5_s4';
eitfile = 't5_s4_inj1';
def = 'Q_R';
% spectfile = 't2_s4';
% eitfile = 't2_s4_inj1';
% def = 'V_L';

do_v = 1;
do_q = 1;
do_vq = 1;

bg = 1; % Background -> 0: Dark, 1: Gray, 2: Light
cl = 1; % Color -> 0: White, 1: Color-White, 2: Color
rg = 1; % Range -> 0: Narrow, 1: Wide
ct = 1; % Contour -> 0: Soft, 1: Hard
md = 2; % Mid -> 0: None, 1: Green, 2: White

bg_cl = [0 0.75 1];
bg = bg_cl(bg+1);
c_bg = bg * [1 1 1]; 

if md == 1
    c_md = 'g';
else 
    c_md = 'w';
end

if rg == 0
    if md == 0
        cm_vq  = buildcmap('krkykckbk') ...
        + 0.5  * buildcmap('rkykckckb') ...
        + 0.5  * buildcmap('kkrkykbkk');
    else
        cm_vq  = buildcmap(['krky',c_md,'ckbk']) ...
        + 0.5  * buildcmap( 'rkykkkckb') ...
        + 0.5  * buildcmap( 'kkrkkkbkk');
    end
   
else
    if md == 0
        cm_vq  = buildcmap('krykkkcbk') ...
        + 0.75 * buildcmap('kkkykckkk') ...    
        + 0.5  * buildcmap('rkkkykkkb') ...
        + 0.5  * buildcmap('kkkkckkkk') ...
        + 0.25 * buildcmap('kkkckykkk');                  

    %         cm_vq   = buildcmap('krkkykkkkkckkbk') ...
%         + 0.875 * buildcmap('kkkkkykkkckkkkk') ...  
%         + 0.75  * buildcmap('kkkykkykckkckkk') ...
%         + 0.5   * buildcmap('rkrkkkkykkkkbkb') ...
%         + 0.5   * buildcmap('kkykkkkckkkkckk') ...
%         + 0.25  * buildcmap('kkkrkkckykkbkkk') ...
%         + 0.125 * buildcmap('kkkkkckkkykkkkk');        
        
    else
        cm_vq  = buildcmap(['kkrkykkk',c_md,'kkkckbkk']) ...
        + 0.75 * buildcmap(['krkkkyk' ,c_md,'k',c_md,'kckkkbk']) ...            
        + 0.5  * buildcmap(['rkkrkk'  ,c_md,'kkk',c_md,'kkbkkb']) ...
        + 0.5  * buildcmap( 'kkkykkykkkckkckkk') ...
        + 0.25 * buildcmap(['kkkkk'   ,c_md,'kykck',c_md,'kkkkk']);
        
%         cm_vq  = buildcmap(['krkkykkk',c_md,          'kkkckkbk']) ...
%         + 0.75 * buildcmap(['kkkykyk' ,c_md,'k',c_md,  'kckckkk']) ...    
%         + 0.5  * buildcmap(['rkrkkk'  ,c_md,'kkk',c_md, 'kkkbkb']) ...
%         + 0.5  * buildcmap( 'kkykkkykkkckkkckk') ...
%         + 0.25 * buildcmap(['kkkrk'   ,c_md,'kykck',c_md,'kbkkk']);     
    
    end
end

cm_vq = min(cm_vq,1);

if cl == 0 % White
    if ct == 0
        cm_v = buildcmap('kbw') + bg * buildcmap('wkk');
        cm_q = buildcmap('krw') + bg * buildcmap('wkk');
    else
        cm_v = buildcmap('kbw') + 0.5 * buildcmap('bkk');
        cm_q = buildcmap('krw') + 0.5 * buildcmap('rkk');
    end
    
elseif cl == 1 % Color-White
    if ct == 0
        cm_v = buildcmap('kbcw') + bg * buildcmap('wkkk');
        cm_q = buildcmap('kryw') + bg * buildcmap('wkkk');
    else
        cm_v = buildcmap('kbcw') + 0.5 * buildcmap('bkkk');
        cm_q = buildcmap('kryw') + 0.5 * buildcmap('rkkk');
    end
    
elseif cl == 2 % Color
    if ct == 0
        cm_v = buildcmap('kbc') + bg * buildcmap('wkk');
        cm_q = buildcmap('kry') + bg * buildcmap('wkk');
    else
        cm_v = buildcmap('kbc') + 0.5 * buildcmap('bkk');
        cm_q = buildcmap('kry') + 0.5 * buildcmap('rkk');
    end
    
else
    error('Invalid color.');
    
end

%%

% V/Q
if do_vq
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\spect\img\spect_vq\low\def_spect_img-',spectfile,'-spect_vq-',def,'-low.fig'],1);
    colormap(cm_vq);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img\eit_vq\spect\def_eit_img-',eitfile,'-eit_vq-',def,'.fig'],1)
    colormap(cm_vq);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img_corr\eit_vq\spect\def_eit_img_corr-',eitfile,'-eit_vq-',def,'.fig'],1)
    colormap(cm_vq);
    set(gca,'Color',c_bg);
end

% V
if do_v
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\spect\img\spect_v\low\def_spect_img-',spectfile,'-spect_v-',def,'-low.fig'],1);
    colormap(cm_v);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img\eit_v\spect\def_eit_img-',eitfile,'-eit_v-',def,'.fig'],1)
    colormap(cm_v);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img_corr\eit_v\spect\def_eit_img_corr-',eitfile,'-eit_v-',def,'.fig'],1)
    colormap(cm_v);
    set(gca,'Color',c_bg);
end

% Q
if do_q
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\spect\img\spect_q\low\def_spect_img-',spectfile,'-spect_q-',def,'-low.fig'],1);
    colormap(cm_q);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img\eit_q\spect\def_eit_img-',eitfile,'-eit_q-',def,'.fig'],1)
    colormap(cm_q);
    set(gca,'Color',c_bg);
    uiopen(['D:\Versuche\Uppsala\Auswertung\work\data\def\eit\img_corr\eit_q\spect\def_eit_img_corr-',eitfile,'-eit_q-',def,'.fig'],1)
    colormap(cm_q);
    set(gca,'Color',c_bg);
end