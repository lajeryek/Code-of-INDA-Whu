function y = qfunc_local(x)
% 本地 Q 函数封装：优先使用内置 qfunc，没有则用 erfc 代替
if exist('qfunc', 'file') == 2
    y = qfunc(x);
else
    y = 0.5 * erfc(x ./ sqrt(2));
end
end