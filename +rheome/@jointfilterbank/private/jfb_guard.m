function jfb_guard(obj, dims, what)
% JFB_GUARD  Refuse an allocation larger than MaxBytes, saying how large it would be.
%
%   jfb_guard(obj, [K nOmega Nf], 'the whole bank')
%
% The point is that an 8.4 GB request fails HERE, at the call site, with a number in the
% message -- not later in the swapper.
%
% Author: Diellor Basha, 2026

    bytes = prod(double(dims)) * 8;
    if bytes > obj.MaxBytes
        error('jointfilterbank:tooLarge', ...
            ['Materialising %s would need %s (%s doubles), over MaxBytes = %s. The bank ' ...
             'stores factors precisely so this is never required: use jointfilters(jfb,m) ' ...
             'for one member, or a reduction (framebounds, scalogram) which loops one ' ...
             'member at a time. Pass ''Force'',true or raise MaxBytes to insist.'], ...
            what, jfb_bytes(bytes), mat2str(dims), jfb_bytes(obj.MaxBytes));
    end
end

% Author: Diellor Basha, 2026
