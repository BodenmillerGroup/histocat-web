if (typeof window.URL.createObjectURL === "undefined") {
  window.URL.createObjectURL = () => {
    // Do nothing
    // Mock this function for mapbox-gl to work
  };
}

HTMLCanvasElement.prototype.getContext = () => {
  // return whatever getContext has to return
};
