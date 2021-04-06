import React from "react";
import TitleInfo from "vitessce/components/TitleInfo";
import { ChannelsSettingsView } from "./ChannelsSettingsView";
import { CoordinationScopes } from "vitessce/types";

type ChannelsSettingsViewSubscriberProps = {
  coordinationScopes: CoordinationScopes;
  removeGridComponent(event: any): void;
};

export function ChannelsSettingsViewSubscriber(props: ChannelsSettingsViewSubscriberProps) {
  const { removeGridComponent } = props;
  return (
    <TitleInfo
      title="Settings"
      removeGridComponent={removeGridComponent}
      isScroll={true}
      theme="dark"
      isReady={true}
    >
      <ChannelsSettingsView />
    </TitleInfo>
  );
}
