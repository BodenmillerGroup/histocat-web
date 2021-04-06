import ScatterplotSubscriber from "../components/scatterplot/ScatterplotSubscriber";
import SpatialSubscriber from "../components/spatial/SpatialSubscriber";
import { ChannelsViewSubscriber } from "../../views/groups/projects/viewers/channels/ChannelsViewSubscriber";
import { SlidesViewSubscriber } from "../../views/groups/projects/viewers/slides/SlidesViewSubscriber";
import { ChannelsSettingsViewSubscriber } from "../../views/groups/projects/viewers/settings/ChannelsSettingsViewSubscriber";
import { BlendViewSubscriber } from "../../views/groups/projects/viewers/image/BlendViewSubscriber";

const registry = {
  slides: SlidesViewSubscriber,
  image: BlendViewSubscriber,
  channels: ChannelsViewSubscriber,
  settings: ChannelsSettingsViewSubscriber,
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
