<template>
  <v-card tile class="ma-6 pa-1">
    <v-card-title>
      <v-row no-gutters>
        <v-col>
          <h5 class="text-h5">{{ group.name }}</h5>
          <span class="text-caption"><v-icon small>mdi-link</v-icon> {{ group.url }}</span>
        </v-col>
      </v-row>
    </v-card-title>
    <v-card-text v-if="group.description">
      {{ group.description }}
    </v-card-text>
    <v-card-actions>
      <v-btn
        v-if="isMember || user.is_admin"
        color="primary"
        :to="{ name: 'main-group', params: { groupId: group.id } }"
      >
        Open
      </v-btn>
      <v-btn
        color="primary"
        :to="{ name: 'main-groups-edit', params: { groupId: group.id } }"
      >
        Edit
      </v-btn>
      <v-btn v-if="group.is_open" color="primary" @click="joinGroup()">
        Request Access
      </v-btn>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { IUserProfile } from "@/modules/user/models";
import { IGroup } from "@/modules/group/models";

@Component
export default class GroupCard extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  @Prop({
    type: Object,
    required: true,
  })
  readonly user!: IUserProfile;

  @Prop({
    type: Object,
    required: true,
  })
  readonly group!: IGroup;

  get userIds() {
    return (this.group as any).members
      ? (this.group as any).members.filter((item) => item.is_active).map((member) => member.user_id)
      : [];
  }

  get isMember() {
    return this.userIds.includes(this.user.id);
  }

  async joinGroup() {
    await this.groupContext.actions.joinGroup(this.group.id);
  }
}
</script>
