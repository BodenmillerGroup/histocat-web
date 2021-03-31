import ScatterplotSubscriber from "../components/scatterplot/ScatterplotSubscriber";
import SpatialSubscriber from "../components/spatial/SpatialSubscriber";
import { SlidesView } from "../../views/groups/projects/viewers/slides/SlidesView";
import { BlendView } from "../../views/groups/projects/viewers/image/BlendView";
import { ChannelsView } from "../../views/groups/projects/viewers/channels/ChannelsView";
import { ChannelsSettingsView } from "../../views/groups/projects/viewers/settings/ChannelsSettingsView";

const registry = {
  slides: SlidesView,
  image: BlendView,
  channels: ChannelsView,
  settings: ChannelsSettingsView,
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
