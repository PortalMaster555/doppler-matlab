function playAudio(amps, Fs, dur)
%PLAYAUDIO Plays the audio vector at sample rate Fs for appx. duration dur
% Exact timing is difficult due to lag in audio player, but it is not 
% needed in this program.
    fprintf("Playing sample of audio...\n")
    audio = audioplayer(amps, Fs); 
    play(audio); 
    pause(dur); 
    stop(audio);
    fprintf("Playback complete.\n")
end

