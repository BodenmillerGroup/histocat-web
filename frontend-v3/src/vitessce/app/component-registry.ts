import ScatterplotSubscriber from '../components/scatterplot/ScatterplotSubscriber';
import SpatialSubscriber from '../components/spatial/SpatialSubscriber';

const registry = {
  scatterplot: ScatterplotSubscriber,
  spatial: SpatialSubscriber,
} as any;

export function getComponent(name: string) {
  const component = registry[name];
  if (component === undefined) {
    throw new Error(`Could not find definition for "${name}" in registry.`);
  }
  return registry[name];
}
