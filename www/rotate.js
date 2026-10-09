// Client-side "Rotate" animation for a slider that drives a server plot.
//
// The next angle is requested only after the previous frame of the plot has
// arrived, so there is never more than one frame in flight and Pause takes
// effect immediately (a server-side timer queues frames faster than ggplot can
// draw them, which made Pause lag by several seconds).
function svRotate(opts) {
  var step = opts.step || 2, max = opts.max || 180;
  var playing = false, waiting = false;

  function $btn() { return $('#' + opts.button); }

  function setPlaying(on) {
    playing = on;
    $btn().html(on ? '<i class="fa fa-pause"></i> Pause' : '<i class="fa fa-play"></i> Rotate');
    if (on) advance();
  }

  function advance() {
    if (!playing || waiting) return;
    var slider = $('#' + opts.slider).data('ionRangeSlider');
    if (!slider) return;
    var v = slider.result.from + step;
    if (v > max) v = 0;
    slider.update({ from: v });
    waiting = true;
    Shiny.setInputValue(opts.slider, v);
  }

  function frameDone(e) {
    if (e.name !== opts.plot) return;
    waiting = false;
    if (playing) setTimeout(advance, 40);
  }

  $(document).on('click', '#' + opts.button, function () {
    setPlaying(!playing);
  });
  (opts.stop || []).forEach(function (id) {
    $(document).on('click', '#' + id, function () { setPlaying(false); });
  });
  $(document).on('shiny:value shiny:error', frameDone);
  // Plot is hidden (other tab) -> it won't re-render, so stop.
  $(document).on('shown.bs.tab', function () { if (playing) setPlaying(false); });
}
