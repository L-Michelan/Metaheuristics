% Função alvo
function y = f_vet(x)
    y = arrayfun(@(x) 2*exp(x) - x*sin(x+3) - 3, x);
end
