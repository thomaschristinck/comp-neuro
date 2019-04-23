function y = two_param_sigmoid(x, params)
% Sigmoidal activation function (controlled by two parameters)
    b = params(1);
	x0 = params(2);
	
	y = 1 ./ (1 + exp(-b .* (x - x0)));
end