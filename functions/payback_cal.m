function [year,npv] = payback_cal(infra,reduction,discount_factor,max_payback,plot_option)
if infra <= 0
    year = 0;
    npv = 0;
else
    npv = infra;
    year = max_payback+1;
    payback = 0;
    for t = 1:1:max_payback
        npv = npv - reduction(t)*(1-discount_factor)^t;
        current(t) = npv;
        if npv <= 0 && payback == 0
            year = t;
            payback = 1;
        end
    end
end

if plot_option >= 1
    figure;
    bar(1:1:max_payback,current);
    xlim([0,max_payback+1])
    xlabel('year');
    ylabel('cost ($)');
end