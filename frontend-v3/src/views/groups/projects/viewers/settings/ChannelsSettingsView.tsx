import { useProjectsStore } from "modules/projects";
import { BrushableHistogram } from "components/BrushableHistogram";

type ChannelsSettingsViewProps = {};

export function ChannelsSettingsView(props: ChannelsSettingsViewProps) {
  const selectedChannels = useProjectsStore((state) => state.getSelectedChannels());

  return (
    <div>
      {selectedChannels.map((channel) => (
        <BrushableHistogram key={channel.name} channel={channel} />
      ))}
    </div>
  );
}
