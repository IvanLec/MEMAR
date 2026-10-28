function PR(score_met, score_not_met; ini=0, pass=0.0001,fin=1)
    """ Function to compute the Prescision Recall curve. It takes as input:
        1) score_met: the probability assigned by a model to real methylated sites
        2) score_not_met: the probability assigned by a model to real non-methylated sites
    Some optional arguments can be passed:
        3) ini: initial treshold value (its default value is 0) 
        4) pass: the variation of the treshold at each iteration (its default value is 0.0001)
        5) fin: last treshold value (its default values is 1)
    The function returns:
        1) R: the Recall
        2) P: the Precision
    """
    #Setting the range of the treshold, from the intial to the final value with the given pass
    range=vcat(ini:pass:fin)
    P, R=[],[] #initialization of empty vectors for precision and recall 
    for t in range #looping over treshold values
        #True Positive 
        tp=count(x->x>=t, score_met) #all predicted sites in the methylated ones which are equal or above the treshold value
        #False Poisitve
        fp=count(x->x>=t, score_not_met) #all predicted sites in the non methylated sites which are equal or above the treshold
        #precision=(True Positive)/(True Positive+ False Positive)
        push!(P, tp/(tp+fp))
        #False Negative
        fn=count(x->x<t, score_met) #sites which are predicted negative (below the treshold value) but they are false (since they are in the methylated dataset)
        #Recall=(True Positive)/(True Positive + False Negative)
        push!(R, tp/(tp+fn))
    end 
    return  R, P
end

function compute_F1_threshold(score_met, score_not_met, threshold)
    """Function to compute the F1 score at a given treshold. It takes an inputs:
        1) score_met: the probability assigned by a model to real methylated sites
        2) score_not_met: the probability assigned by a model to real non-methylated sites
        3) treshold: the value at which the F1 is wanted to be compute_F1_threshold
        The function returns: 
        1) xx: the recall at the given treshold
        2) yy: the precision at given treshold 
        3) F1: the F1 score at given treshold
    """
    #Computing the precision recall just over 1 treshold value
    xx, yy=PR(score_met, score_not_met; ini=threshold,  fin=threshold)
    #Saving the first value corresponding to the chosen treshold value for precision and recall
    xx,yy=xx[1], yy[1]
    #F1 score at the given treshold value
    F1=2 * (yy*xx)./(yy+xx)
    return xx, yy, F1
end

function plot_PR_score(score_met, score_not_met, col, leg; pas=0.001, xx=0.05, yy=xx, ls="-", verbose=false, mark="o")
    """Function to plot Precision Recal Curves. It takes as argument:
        1) score_met: the probability assigned by a model to real methylated sites
        2) score_not_met: the probability assigned by a model to real non-methylated sites
        3) col: the coulor of the precision recall curve 
        4) leg: legend for the given precision recall curve
    Some additional parameters can be passed:
        5) pass: the variation of the treshold at each iteration (its default value is 0.001)
        6) xx: offeset of the x axis to plot a text (default value equal to 0.05)
        7) yy: offeset of the y axis to plot a text (default value equivalent to xx)
        8) ls: linestyle (default value is continous curve '-')
        9) verbose: booelan variable to control functions output (default value is 'false')
        10) mark: marker shape (the default value is a circle 'o')
    The function returns: 
        1) the Precision Recall curve plot
        2) tre2: the treshold value corrresponding to the point with maximum F1 score
        3) area: the AUPRC if verbose is 'false'
        4) mF1: maximum F1 one value if verbose is 'false'
    """
    #Saving Precision Recall Curve
    x, y=PR(score_met, score_not_met; pass=pas)
    #Filtering if NaN values are retrieved over extreme treshold values
    ind=length(x)
    for i in 1:length(y)
        if isnan(y[i])==true
            ind=i-1
            break
        end
    end
    #Saving before extreme values (if present)
    #Setting precision recall starting point by convention  (0, 1)
    x=vcat(x[1:ind],0 ) 
    y=vcat(y[1:ind], 1)
    #Computing AUPRC
    #x=x[1:findlast(x->x>0, x)]
    #y=y[1:findlast(x->x>0, x)]
    #Computing AUPRC
    area = abs(sum((y[i] + y[i+1]) * (x[i+1] - x[i]) / 2 for i in 1:length(x)-1))
    #L=round(digits=3,length(score_met)/(length(score_met)+length(score_not_met))*100)
    #Plotting the result
    plot(x,y, label=leg*"\nAUPRC= $(round(digits=3, (area)) )", color=col, linestyle=ls, linewidth=2)
    #Finding the value that maximize the F1 score
    T=vcat(0:pas:1...)
    F1=2 .* (y.*x)./(y.+x)
    F1=vcat([(isnan(x)==true ? 0 : x) for x in F1]...) #not considering the artificial point (0, 1) since the harmonic mean is not definable (denominator must be different than 0)
    pos_F1=argmax(F1)
    mF1=maximum(F1)
    tre2=(T[pos_F1])
    #plotting the max F1 score and its value
    scatter(x[pos_F1], y[pos_F1], c=col, edgecolor="black", linewidth=0.4, zorder=2, marker=mark, s=80)
    text(x[pos_F1]+xx, y[pos_F1]+yy, "Max F1:$(round(digits=3, (mF1)))", fontsize=12, color="black", zorder=2)
    if verbose==false
        return gcf(), tre2
    else
        return gcf(), tre2, area, mF1
    end
end


function ROC(score_met, score_not_met; ini=0,  pass=0.0001, fin=1)
    """ Function to compute the Receiver operating characteristic curve. It takes as input:
        1) score_met: the probability assigned by a model to real methylated sites
        2) score_not_met: the probability assigned by a model to real non-methylated sites
    Some optional arguments can be passed:
        3) ini: initial treshold value (its default value is 0) 
        4) pass: the variation of the treshold at each iteration (its default value is 0.0001)
        5) fin: last treshold value (its default values is 1)
    The function returns:
        1) N: the False Positive Rate
        2) T: the True Positive Rate
    """
    #Setting the range of the treshold, from the intial to the final value with the given pass
    range=vcat(ini:pass:fin)
    T, N=[],[] #initialization of empty vectors for FPR and TPR 
    for t in range 
        #Recall=TPR=(True Positive)/(True Positive + False Negative)
        push!(T, count(x->x>=t, score_met)/length(score_met))
        #FPR=(False Positive)/(False Positive + True Negative)
        push!(N, count(x->x>=t, score_not_met)/length(score_not_met)) 
    end 
    return  N, T
end

function plot_ROC_score(score_met, score_not_met, col, leg; pas=0.001, xx=0.05, yy=xx, ls="-")
    """Function to plot Receiver Operating Charcteristic Curve. It takes as argument:
        1) score_met: the probability assigned by a model to real methylated sites
        2) score_not_met: the probability assigned by a model to real non-methylated sites
        3) col: the coulor of the precision recall curve 
        4) leg: legend for the given precision recall curve
    Some additional parameters can be passed:
        5) pass: the variation of the treshold at each iteration (its default value is 0.001)
        6) xx: offeset of the x axis to plot a text (default value equal to 0.05)
        7) yy: offeset of the y axis to plot a text (default value equivalent to xx)
        8) ls: linestyle (default value is continous curve '-')
    The function returns: 
        1) the ROC plot
    """
    "Computing the ROC curve"
    x, y=ROC(score_met, score_not_met; pass=pas)
    #Filtering if NaN values are retrieved over extreme treshold values
    ind=length(x)
    for i in 1:length(y)
        if isnan(y[i])==true
            ind=i-1
            break
        end
    end
    x=x[1:ind]
    y=y[1:ind]
    #COmputing AUROC
    area = abs(sum((y[i] + y[i+1]) * (x[i+1] - x[i]) / 2 for i in 1:length(x)-1))
    #L=round(digits=3,length(score_met)/(length(score_met)+length(score_not_met))*100)
    #Plotting ROC
    plot(x,y, label=leg*"\nAUROC= $(round(digits=3, (area)) )", color=col, linestyle=ls, linewidth=2)
    return gcf()
end
