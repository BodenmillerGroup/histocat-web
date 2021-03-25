export const getCursorWithTool = () => 'crosshair';
export const getCursor = (interactionState: { isDragging: any; }) => (interactionState.isDragging
  ? 'grabbing' : 'default'
);
